library(ggplot2)
library(dplyr)

phi <- function(n, m, j, k)
{
  1 - (j+n*k) / (m * (n^2+n))  
}

grid <- expand.grid(n=1:30,m=1:5,j=1,k=1:3)
grid$phi = sapply(1:nrow(grid), function(i) {
  phi(n=grid$n[i],
      m=grid$m[i],
      j=grid$j[i],
      k=grid$k[i])
  })

g = ggplot(grid %>%
         mutate(m=factor(m,ordered=TRUE), k=factor(k,ordered=TRUE)), 
       aes(x=n,y=phi,color=m,group=paste(m,k),linetype=k)) + 
  geom_line() +
  theme_bw() +
  xlab(expression(paste("Number of species (", italic(n),")"))) +
  ylab(expression(paste("Fractional parameter reduction (", italic(phi),")"))) +
  geom_hline(color='red',yintercept = 0) +
  annotate(geom='text',x=30,y=-1,label=expression(paste(italic(j),"=1 environment variables")),adj=1,size=3) +
  scale_color_viridis_d(name=expression(paste("Number of environments (", italic(m), ")"))) +
  scale_linetype_discrete(name=expression(paste("Number of traits (", italic(k), ")"))) +
  theme(legend.position = 'right')

ggsave(g, file='parameter_reduction.pdf',width=6,height=4)

       