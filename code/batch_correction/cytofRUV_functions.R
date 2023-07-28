# CytofRUV Normalise data functions.
# These are not exported out of the package, so we can't just call them
# directly from the package.
# Why not just use normalise_data?
# Because it requires untransformed FCS files stored as flowFrames (flowSet)
# which is very hard (maybe impossible to do) if we use supercell since
# we do transformation then run supercell.


### Functions
make_residual_mat_several_rep<- function(data,clusters, norm_clus_list,samples,rep_samples_list){
    #make_residual_mat(raw_Y,data$cluster, norm_clusters,norm_clusters_second,data$sample,rep_samples,second_rep_samples)
    res_mat=data
    res_mat[]=0
    for (r in 1:length(rep_samples_list)){
        norm_clus=norm_clus_list[[r]]
        rep_samples=rep_samples_list[[r]]
        mean_pseud=matrix(nrow = length(norm_clus),ncol=dim(data)[2])
        for (i in 1:length(norm_clus)){
            tmp=((clusters == norm_clus[i])&(samples%in%rep_samples))
            mean_pseud[i,]=colMeans(data[tmp,])
            res_mat[tmp,]<- t(apply(data[tmp,], 1, function(x) x-mean_pseud[i,]))
        }
    }
    return(res_mat)
}


fastRUVIII = function(Y, M, ctl,res_mat,k=NULL, eta=NULL, average=FALSE, fullalpha=NULL){
    # Assumes good input
    if (!(k > 0)) stop("Bad input - read the documentation")
    Y = ruv::RUV1(Y,eta,ctl)
    m = nrow(Y)
    #Y0 = fast_residop(Y, M)
    Y0=res_mat
    fullalpha = diag(rsvd(Y0)$d) %*% t(rsvd(Y0)$v)
    alpha = fullalpha[1:k,,drop=FALSE]
    ac = alpha[,ctl,drop=FALSE]
    W = Y[,ctl] %*% t(ac) %*% solve(ac %*% t(ac))
    newY = Y - W %*% alpha
    return(list(newY = newY, fullalpha=fullalpha))
}


fast_residop <- function (A, B) {
    return(A - B %*% solve(t(B) %*% B) %*% (t(B) %*% A))
}


run_RUVIII <- function(data, norm_clusters, k,rep_samples){
    raw_Y <- as.matrix(data[3:ncol(data)])
    # Standardise the input and then compensate output
    col_means <- colMeans(raw_Y)
    col_sds <- apply(raw_Y, 2, function(x) sd(x))

    for(i in 1:ncol(raw_Y)){
        raw_Y[,i] <- (raw_Y[,i] - col_means[i])/col_sds[i]
    }
    # Run the actual RUVIII
    res_mat<-make_residual_mat_several_rep(raw_Y,data$cluster, norm_clusters,data$sample,rep_samples)
    norm_Y <- fastRUVIII(Y = raw_Y, M, ctl = c(1:ncol(raw_Y)), res_mat=res_mat,k = k)$newY

    for(i in 1:ncol(norm_Y)){
        norm_Y[,i] <- norm_Y[,i]*col_sds[i] + col_means[i]
    }

    return(norm_Y)
}
