########################## Just a network plot
set.seed(46+2)

V = 1            # One blocking variable
G = 3            # Three categories in this variable
N_id = 65        # Number of people

clique = sample(1:3, N_id, replace=TRUE)
B = matrix(-10, nrow=G, ncol=G) 
diag(B) = -5.28 # Block matrix

B[1,3] = -10.9
B[3,2] = -9.9

B = B-1

A = simulate_sbm_plus_srm_network(N_id = N_id, B=list(B=B), V=V, 
                         groups=data.frame(clique=factor(clique)),
                         individual_predictor=matrix(rnorm(N_id, 0, 1), nrow=N_id, ncol=1), 
                         individual_effects=matrix(c(0.7, 0.9),ncol=1, nrow=2),
                         sr_sigma = c(1.8, 1.3), sr_rho = 0.75,
                         dr_sigma = 3.7, dr_rho = 0.97,
                         mode="bernoulli"
                               )

cols = plvs_vltra("mystic_mausoleum")

Net = graph_from_adjacency_matrix(A$network, mode = c("directed"))
V(Net)$color = c(cols[5], cols[9], cols[7])[A$group_ids$clique]

E(Net)$color = c("grey60","black")[is.mutual(Net)+1]

pdf("Example.pdf", height=8, width=8)
par(bg="white")
plot(Net, edge.arrow.size =0.3, edge.curved = 0.3, vertex.label=NA, vertex.size = 5)
dev.off()



