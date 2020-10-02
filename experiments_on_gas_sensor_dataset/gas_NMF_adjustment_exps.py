import os
import numpy as np

# set rwnmf default options
def get_options(opt):
    opt_fields = list(opt.keys())

    if not ('random_state' in opt_fields):
        opt['random_state'] = 0

    if not ('max_iters' in opt_fields):
        opt['max_iters'] = 500

    opt['tol_fit_improvement'] = 1e-4
    opt['tol_fit_error'] = 1e-4

    return opt


# regularized non-negative matrix factorization
def rwnmf(X, k, alpha, options=None):
    # applies regularized weighted nmf to matrix X with k factors
    # ||X-UV^T||
    # set default options if needed
    if options is None:
        options = {}
    options = get_options(options)
    eps = np.finfo(np.float).eps
    early_stop = False

    # get observations matrix W
    W = np.isnan(X)
    X[W] = 0  # set missing entries as 0
    W = ~W

    # initialize factor matrices
    rnd = np.random.RandomState(options['random_state'])
    U = rnd.rand(X.shape[0], k)
    U = np.maximum(U, eps)

    V = np.linalg.lstsq(U, X, rcond=None)[0].T
    V = np.maximum(V, eps)

    Xr = np.inf * np.ones(X.shape)

    for i in range(options['max_iters']):
        # update U
        U = U * np.divide(((W * X) @ V), (W * (U @ V.T) @ V + alpha * U))
        U = np.maximum(U, eps)
        # update V
        V = V * np.divide((np.transpose(W * X) @ U), (np.transpose(W * (U @ V.T)) @ U + alpha * V))
        V = np.maximum(V, eps)

        # compute the resduals
        if i % 10 == 0:
            # compute error of current approximation and improvement
            Xi = U @ V.T
            fit_error = np.linalg.norm(X - Xi, 'fro')
            fit_improvement = np.linalg.norm(Xi - Xr, 'fro')

            # update reconstruction
            Xr = np.copy(Xi)

            # check if early stop criteria is met
            if fit_error < options['tol_fit_error'] or fit_improvement < options['tol_fit_improvement']:
                error = fit_error
                early_stop = True
                break

    if not early_stop:
        Xr = U @ V.T
        error = np.linalg.norm(X - Xr, 'fro')

    return Xr, U, V, error


def get_linear_extension(X, sequence,Xx):
    # NOTE: this function should only be executed if there are no observed values for any house in the gap start/end
    Yr = np.copy(X)

    non_obs_idxs=list(np.sum(Xx[:,sequence]!=Xx[:,sequence], 0)==Xx.shape[0])
    obs_sequence = [i for i in range(len(sequence)) if not(non_obs_idxs[i])]
    non_obs_sequence = [i for i in range(len(sequence)) if non_obs_idxs[i]]

    # if both borders exist consider linear interpolation
    x_left = Yr[obs_sequence[0]]
    x_right = Yr[obs_sequence[-1]]
    delta_x = (x_right - x_left) / (len(obs_sequence)-1)

    # extend left linearly
    if non_obs_sequence[0]==0:
        left_idx=[non_obs_sequence[0]]
        ctr = 1 # index to second non observed position
        while ctr<len(non_obs_sequence) and non_obs_sequence[ctr]==non_obs_sequence[ctr-1]+1:
            left_idx.append(non_obs_sequence[ctr])
            ctr+=1

        for i in reversed(range(len(left_idx))):
            Yr[left_idx[i]] = x_left + (i-len(left_idx)) * delta_x

    # extend right linearly
    if non_obs_sequence[-1]==len(Yr)-1:
        right_idx=[non_obs_sequence[-1]]
        ctr = len(non_obs_sequence)-2 # index to penultimate non observed position
        while ctr>=0 and non_obs_sequence[ctr]==non_obs_sequence[ctr+1]+-1:
            right_idx.append(non_obs_sequence[ctr])
            ctr-=1
        right_idx.reverse()
        for i in range(len(right_idx)):
            Yr[right_idx[i]] = x_left + (len(obs_sequence)+i) * delta_x

    return Yr


def get_linear_reco(X, interval, sequence):
    Yr = np.copy(X)
    x_left = Yr[interval[0], sequence[0] - 1]
    x_right = Yr[interval[0], sequence[-1] + 1]

    delta_x = (x_right - x_left) / (len(sequence) + 1)
    for i in range(len(sequence)):
        Yr[interval[0], sequence[i]] = x_left + (i + 1) * delta_x

    return Yr[interval[0], sequence]


def get_nmf_reco(X, k, alpha, options):
    Y = np.copy(X)

    # scale data
    Y, scaling_params = scale_data(Y)

    # apply rNMF
    out = rwnmf(Y, k, alpha, options)
    Xr = out[0]

    # re-scale output
    for i in range(X.shape[0]):
        Xr[i,:] = Xr[i,:] * (scaling_params[i][1] - scaling_params[i][0]) + scaling_params[i][0]

    return Xr


def get_aligned_nmf_reco(X, sequence, Xr, Yr):
    # extend reconstruction with linear trend
    if all(X[:, sequence[0]]!=X[:, sequence[0]]) or all(X[:, sequence[-1]]!=X[:, sequence[-1]]):
        Xr = get_linear_extension(Xr, sequence, X)
    Xr_shifted = np.copy(Xr)
    
    # get "linear" trend
    Zr = np.zeros((len(sequence)))
    x_left = Xr_shifted[0]
    x_right = Xr_shifted[-1]
    delta_x = (x_right - x_left) / (len(sequence) + 1)
    for i in range(len(sequence)):
        Zr[i] = x_left + (i + 1) * delta_x

    # align
    Xr_shifted = Xr_shifted - (Zr - Yr)

    return Xr_shifted


def load_data():
    X=np.load('data/gas_sensor.npy') 
    
    return X


def scale_data(X):
    Y=np.copy(X)
    scaling_params = []
    for i in range(Y.shape[0]):
        y=Y[i,Y[i,]==Y[i,]]
        if len(y)>0 and min(y)!=max(y):
            scaling_params.append([min(y),max(y)])
            Y[i,]=(Y[i,]-min(y))/(max(y)-min(y))
        else:
            scaling_params.append([0,0])

    return Y, scaling_params


def corrupt_data(X, gap_size, nblocks, sample,seed):
    Y = np.copy(X)

    np.random.seed(seed+len(sample))

    # create sample
    for nblock in range(nblocks):
        i = np.random.choice(X.shape[0], 1)[0]
        j = np.random.randint(0, X.shape[1] - gap_size - 1, 1)[0]
        # avoid overlaps
        while any([j in list(range(interval[1]-1,interval[1]+gap_size+1)) for interval in sample if interval[0]==i]) or any([j+gap_size-1 in list(range(interval[1]-1,interval[1]+gap_size+1)) for interval in sample if interval[0]==i]):
            i = np.random.choice(X.shape[0], 1)[0]
            j = np.random.randint(0, X.shape[1] - gap_size - 1, 1)[0]

        Y[i, list(range(j,j+gap_size))]=np.nan

        # store details
        sample.append([i, j])

    return Y, sample


def sequence_error(x, y):
    error = 0
    for i in range(len(x)):
        error += np.power(x[i] - y[i], 2)

    return error


def main():
    random_state=3445
    alpha = 0.1
    Kmax=12
    options = {}
    options['random_state'] = random_state
    missing_rate=0.5
    sampleseeds = [71049, 33, 8101, 59174, 223344, 1111, 431, 24560, 6712, 9999]
    gap_sizes = [5,10,15,20,30]
    dataset ='gas_sensor'

    print('########################################')
    print('EXPERIMENTS ON '+dataset.upper())
    print('########################################')

    # load dataset
    X=load_data()
    
    # initialize results struct
    press_nmf={}
    press_adjusted_nmf={}
    for gap_size in gap_sizes:
        press_nmf[gap_size] = {}
        press_adjusted_nmf[gap_size] = {}
        press_nmf[gap_size][missing_rate]=np.zeros((len(sampleseeds),Kmax))
        press_adjusted_nmf[gap_size][missing_rate] = np.zeros((len(sampleseeds), Kmax))

    ctr = 0
    for sampleseed in sampleseeds:
        for gap_size in gap_sizes:
            Y=np.copy(X)

            # set the number of blocks according to the number of timestamps
            nblocks = int(np.round(X.shape[1]/gap_size*X.shape[0])) # total number of gaps with no overlap
            nblocks=int(np.floor(nblocks*missing_rate))

            sample= []
            # corrupt data
            Y, sample =corrupt_data(Y,gap_size, nblocks,sample,sampleseed)

            # apply cp with increasing number of components
            for k in range(1, Kmax + 1):
                # get nmf reconstruction
                Xr = get_nmf_reco(Y, k, alpha, options)

                # adjust each of the gaps
                for interval in sample:
                    # set indexes in the season context
                    sequence = list(range(interval[1], interval[1] + gap_size))

                    # get cp reconstruction
                    Yr = get_linear_reco(Y, interval, sequence)

                    # get press
                    nmf_error = sequence_error(X[interval[0], sequence], Xr[interval[0], sequence])
                    press_nmf[gap_size][missing_rate][ctr, k - 1] += nmf_error

                    # get alignment curve
                    Xr_shifted = get_aligned_nmf_reco(Y, sequence, Xr[interval[0], sequence], Yr)

                    # get error
                    adjusted_nmf_error = sequence_error(X[interval[0], sequence], Xr_shifted)
                    press_adjusted_nmf[gap_size][missing_rate][ctr,k - 1] += adjusted_nmf_error
        ctr+=1

    # get  results summary
    press_nmf_summary=np.zeros(len(gap_sizes))
    press_adjusted_nmf_summary = np.zeros(len(gap_sizes))

    for gap_size in gap_sizes:
        print('----------------------------')
        print('| Gap size of '+str(gap_size)+' timestamps |')
        gap_size_idx=gap_sizes.index(gap_size)
        print('>> Missing rate of ' + str(missing_rate*100) + '%')

        mean_press_nmf = np.mean(press_nmf[gap_size][missing_rate],0)
        mean_adjusted_press_nmf = np.mean(press_adjusted_nmf[gap_size][missing_rate], 0)

        min_press_nmf = np.min(mean_press_nmf)
        min_press_idx = np.argmin(mean_press_nmf )
        min_press_adjusted_nmf = mean_adjusted_press_nmf[min_press_idx]

        press_nmf_summary[gap_size_idx]=min_press_nmf
        press_adjusted_nmf_summary[gap_size_idx] = min_press_adjusted_nmf

        # print results
        print(' - NMF ' + str(min_press_nmf))
        print(' - adjusted NMF ' + str(min_press_adjusted_nmf))

    # store results
    np.savez('results/'+dataset+'_missrate'+str(missing_rate), press_nmf=press_nmf,press_adjusted_nmf=press_adjusted_nmf, press_nmf_summary=press_nmf_summary,press_adjusted_nmf_summary=press_adjusted_nmf_summary)


if __name__ == '__main__':
    main()





