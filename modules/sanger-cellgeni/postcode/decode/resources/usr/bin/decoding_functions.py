import numpy as np
import pyro
import torch
from pyro import poutine
from pyro.distributions import Categorical, constraints, MultivariateNormal
from pyro.infer import config_enumerate, SVI, TraceEnum_ELBO
from pyro.infer.autoguide import AutoDelta
from pyro.optim import Adam
from tqdm import tqdm

assert pyro.__version__.startswith("1")

# import pysnooper


# auxiliary functions required for decoding
def torch_format(numpy_array):
    D = numpy_array.shape[1] * numpy_array.shape[2]
    return (
        torch.tensor(numpy_array)
        .float()
        .transpose(1, 2)
        .reshape(numpy_array.shape[0], D)
    )


def barcodes_01_from_channels(barcodes_1234, C, R):
    K = barcodes_1234.shape[0]
    barcodes_01 = np.ones((K, C, R))
    for b in range(K):
        barcodes_01[b, :, :] = 1 * np.transpose(
            barcodes_1234[b, :].reshape(R, 1) == np.arange(1, C + 1)
        )
    return barcodes_01


def kronecker_product(tr, tc):
    tr_height, tr_width = tr.size()
    tc_height, tc_width = tc.size()
    out_height = tr_height * tc_height
    out_width = tr_width * tc_width
    tiled_tc = tc.repeat(tr_height, tr_width)
    expanded_tr = (
        tr.unsqueeze(2)
        .unsqueeze(3)
        .repeat(1, tc_height, tc_width, 1)
        .view(out_height, out_width)
    )
    return expanded_tr * tiled_tc


def chol_sigma_from_vec(sigma_vec, D):
    L = torch.zeros(D, D)
    L[torch.tril(torch.ones(D, D)) == 1] = sigma_vec
    return torch.mm(L, torch.t(L))


@config_enumerate
def model_constrained_tensor(data, N, D, C, R, K, codes, batch_size=None):
    w = pyro.param("weights", torch.ones(K) / K, constraint=constraints.simplex)

    # using tensor sigma
    sigma_ch_v = pyro.param("sigma_ch_v", torch.eye(C)[np.tril_indices(C, 0)])
    sigma_ch = chol_sigma_from_vec(sigma_ch_v, C)
    sigma_ro_v = pyro.param("sigma_ro_v", torch.eye(D)[np.tril_indices(R, 0)])
    sigma_ro = chol_sigma_from_vec(sigma_ro_v, R)
    sigma = kronecker_product(sigma_ro, sigma_ch)

    # codes_tr_v = pyro.param('codes_tr_v', 3 * torch.ones(1, D), constraint=constraints.positive)
    codes_tr_v = pyro.param(
        "codes_tr_v", 3 * torch.ones(1, D), constraint=constraints.greater_than(1.0)
    )
    codes_tr_consts_v = pyro.param("codes_tr_consts_v", -1 * torch.ones(1, D))

    theta = torch.matmul(codes * codes_tr_v + codes_tr_consts_v, mat_sqrt(sigma, D))

    with pyro.plate("data", N, batch_size) as batch:
        z = pyro.sample("z", Categorical(w))
        pyro.sample("obs", MultivariateNormal(theta[z], sigma), obs=data[batch])


auto_guide = AutoDelta(
    poutine.block(
        model_constrained_tensor,
        expose=[
            "weights",
            "codes_tr_v",
            "codes_tr_consts_v",
            "sigma_ch_v",
            "sigma_ro_v",
        ],
    )
)


def mat_sqrt(A, D):
    try:
        U, S, V = torch.svd(A + 1e-3 * A.mean() * torch.rand(D, D))
    except Exception:
        U, S, V = torch.svd(A + 1e-2 * A.mean() * torch.rand(D, D))
    S_sqrt = torch.sqrt(S)
    return torch.mm(torch.mm(U, torch.diag(S_sqrt)), V.t())


def train(
    svi, num_iterations, data, N, D, C, R, K, codes, print_training_progress, batch_size
):
    pyro.clear_param_store()
    losses = []
    for j in tqdm(range(num_iterations), disable=not print_training_progress):
        loss = svi.step(data, N, D, C, R, K, codes, batch_size)
        losses.append(loss)
    return losses


# input - output decoding function
# @pysnooper.snoop()
def train_and_save_model(
    spots,
    barcodes_01,
    estimate_bkg=True,
    estimate_additional_barcodes=None,
    num_iter=60,
    batch_size=15000,
    up_prc_to_remove: float = 99.95,
    print_training_progress=True,
    set_seed=1,
):
    # INPUT:
    # spots: a numpy array of dim N x C x R (number of spots x coding channels x rounds);
    # barcodes_01: a numpy array of dim K x C x R (number of barcodes x coding channels x rounds)
    # OUTPUT:
    # 'class_probs': posterior probabilities computed via e-step
    # 'class_ind': indices of different barcode classes (genes / background / infeasible / nan)
    # 'params': estimated model parameters
    # 'norm_const': constants used for normalization of spots prior to model fitting

    # if cuda available, runs on gpu
    torch.set_default_tensor_type("torch.FloatTensor")

    N = spots.shape[0]
    if N == 0:
        print("There are no spots to decode.")
        return
    C = spots.shape[1]
    R = spots.shape[2]
    K = barcodes_01.shape[0]
    D = C * R
    data = torch_format(spots)
    codes = torch_format(barcodes_01)

    if estimate_bkg:
        bkg_ind = codes.shape[0]
        codes = torch.cat((codes, torch.zeros(1, D)))
    else:
        bkg_ind = np.empty((0,), dtype=np.int32)
    if np.any(estimate_additional_barcodes is not None):
        inf_ind = codes.shape[0] + np.arange(estimate_additional_barcodes.shape[0])
        codes = torch.cat((codes, torch_format(estimate_additional_barcodes)))
    else:
        inf_ind = np.empty((0,), dtype=np.int32)

    # filter out spots that doesn't have appropriate profile
    ind_keep = (
        np.where(
            np.sum(
                data.cpu().numpy()
                < np.percentile(data.cpu().numpy(), up_prc_to_remove, axis=0),
                axis=1,
            )
            == D
        )[0]
        if up_prc_to_remove < 100
        else np.arange(0, N)
    )
    if len(ind_keep) == 0:
        raise ValueError(
            "No spots remain after profile filtering; please check the input data "
            "or increase up_prc_to_remove."
        )
    filtered_data = data[ind_keep, :]
    # Get statistics for normalization after removing outliers
    s = torch.tensor(np.percentile(filtered_data.cpu().numpy(), 60, axis=0))
    max_s = torch.tensor(np.percentile(filtered_data.cpu().numpy(), 99.9, axis=0))
    # print(s, max_s)
    min_s = torch.min(filtered_data, dim=0).values
    log_add = (s**2 - max_s * min_s) / (max_s + min_s - 2 * s)
    log_add = torch.max(
        -torch.min(filtered_data, dim=0).values + 1e-10, other=log_add.float()
    )
    data_log_added = data + log_add
    if 0 in data_log_added:
        data_log_added[data_log_added == 0] = np.finfo(
            float
        ).eps  # Replace values 0 with the smallest positive float
    data_log = torch.log10(data_log_added)
    data_log_mean = data_log[ind_keep, :].mean(dim=0, keepdim=True)
    data_log_std = data_log[ind_keep, :].std(dim=0, keepdim=True)
    data_norm = (data_log - data_log_mean) / data_log_std  # column-wise normalization

    # model training:
    optim = Adam({"lr": 0.085, "betas": [0.85, 0.99]})
    svi = SVI(
        model_constrained_tensor,
        auto_guide,
        optim,
        loss=TraceEnum_ELBO(max_plate_nesting=1),
    )
    pyro.set_rng_seed(set_seed)
    losses = train(
        svi,
        num_iter,
        data_norm[ind_keep, :],
        len(ind_keep),
        D,
        C,
        R,
        codes.shape[0],
        codes,
        print_training_progress,
        min(len(ind_keep), batch_size),
    )
    # collect estimated parameters
    w_star = pyro.param("weights").detach()
    sigma_ch_v_star = pyro.param("sigma_ch_v").detach()
    sigma_ro_v_star = pyro.param("sigma_ro_v").detach()
    sigma_ro_star = chol_sigma_from_vec(sigma_ro_v_star, R)
    sigma_ch_star = chol_sigma_from_vec(sigma_ch_v_star, C)
    sigma_star = kronecker_product(sigma_ro_star, sigma_ch_star)
    codes_tr_v_star = pyro.param("codes_tr_v").detach()
    codes_tr_consts_v_star = pyro.param("codes_tr_consts_v").detach()
    theta_star = torch.matmul(
        codes * codes_tr_v_star + codes_tr_consts_v_star, mat_sqrt(sigma_star, D)
    )

    print(bkg_ind, inf_ind, estimate_bkg)

    return {
        "w_star": w_star,
        "sigma_ch_v_star": sigma_ch_v_star,
        "sigma_ro_v_star": sigma_ro_v_star,
        "sigma_ro_star": sigma_ro_star,
        "sigma_ch_star": sigma_ch_star,
        "sigma_star": sigma_star,
        "codes_tr_v_star": codes_tr_v_star,
        "codes_tr_consts_v_star": codes_tr_consts_v_star,
        "theta_star": theta_star,
        "losses": losses,
        "data_norm": data_norm.cpu(),  # Ensure data_norm is saved as a CPU tensor
        "K": K,
        "C": C,
        "D": D,
        "N": N,
        "R": R,
        "codes": codes.cpu(),  # Ensure codes is saved as a CPU tensor
        "norm_log_add": log_add,
        "norm_data_log_mean": data_log_mean,
        "norm_data_log_std": data_log_std,
        "bkg_ind": np.array(bkg_ind),
        "inf_ind": np.array(inf_ind),
        "estimate_bkg": estimate_bkg,
    }
