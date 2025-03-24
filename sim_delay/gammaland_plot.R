library(ggplot2)
library(tidyverse)
library(patchwork)
rm(list=ls())
set.seed(1)
#First the old point that the same distribution at a different length scale assigns lower likelihoods
df1<- data.frame(x = seq(0,100, length.out = 101)) %>%
  mutate(y1_1 = dgamma(x, shape = 1, scale = 1),
         y1_2 = dgamma(x, shape = 1, scale = 2),
         y1_4 = dgamma(x, shape = 1, scale = 4),
         y1_8 = dgamma(x, shape = 1, scale = 8),
         y1_16 = dgamma(x, shape = 1, scale = 16),
         y1_32 = dgamma(x, shape = 1, scale = 32))

df2<- data.frame(x = seq(0,100, length.out = 101)) %>%
  mutate(y2_1 = dgamma(x, shape = 2, scale = 1),
         y2_2 = dgamma(x, shape = 2, scale = 2),
         y2_4 = dgamma(x, shape = 2, scale = 4),
         y2_8 = dgamma(x, shape = 2, scale = 8),
         y2_16 = dgamma(x, shape = 2, scale = 16),
         y2_32 = dgamma(x, shape = 2, scale = 32))


df3<- data.frame(x = seq(0,100, length.out = 101)) %>%
  mutate(y4_1 = dgamma(x, shape = 4, scale = 1),
         y4_2 = dgamma(x, shape = 4, scale = 2),
         y4_4 = dgamma(x, shape = 4, scale = 4),
         y4_8 = dgamma(x, shape = 4, scale = 8),
         y4_16 = dgamma(x, shape = 4, scale = 16),
         y4_32 = dgamma(x, shape = 4, scale = 32))

df.l1<-df1 %>% gather(key, y, y1_1:y1_32) %>% mutate(key.f = factor(key, levels = c('y1_1', 'y1_2', 'y1_4', 'y1_8', 'y1_16', 'y1_32')))
df.l2<-df2 %>% gather(key, y, y2_1:y2_32) %>% mutate(key.f = factor(key, levels = c('y2_1', 'y2_2', 'y2_4', 'y2_8', 'y2_16', 'y2_32')))
df.l3<-df3 %>% gather(key, y, y4_1:y4_32) %>% mutate(key.f = factor(key, levels = c('y4_1', 'y4_2', 'y4_4', 'y4_8', 'y4_16', 'y4_32')))

ggplot(df.l1, aes(x = x, y = y, colour = key.f)) +
  geom_line()

ggplot(df.l2, aes(x = x, y = y, colour = key.f)) +
  geom_line()

ggplot(df.l3, aes(x = x, y = y, colour = key.f)) +
  geom_line()


#So same distribution under multiplicative transformation = everything becomes less likely

#What about the posterior-predictive, marginalizing over possible values of k and theta given some evidence.... 
#TBC need to finesse / search for the right close-to-uniform prior to avoid NaNs...

#Hyper prior params that I think are fairly close to uninformative
p = 1#0
q = 1#0
# r = 000000.1#0
# s = 000000.1#0
r = 0.000001#0
s = 0.000001#0

x1<-c(20,20,20,20,20)
x2<-c(100,100,100,100,100)

#Posterior params after seeing x1
p_1 = p*prod(x1)
q_1 = q+sum(x1)
r_1 = r+length(x1)
s_1 = s+length(x1)

#Posterior params after seeing x2
p_2 = p*prod(x2)
q_2 = q+sum(x2)
r_2 = r+length(x2)
s_2 = s+length(x2)

#A grid of values of k and theta...
ix<-expand.grid(k = seq(0.01,20, length.out = 100), theta = seq(0.01,20,length.out = 100))

#The prior probability at each grid location approximating p(k,theta)
ix$y.prior.un<-(p^(ix$k-1)*exp(-ix$theta^-1*q)) / (gamma(ix$k)^r * ix$theta^{ix$k*s})
#Normalize it to sum to 1 since we cannot compute Z exactly
ix$y.prior<-ix$y.prior.un/sum(ix$y.prior.un, na.rm=T)

#The posterior probability at each grid location approximating p(k,theta|x1)
ix$y.un1<-(p_1^(ix$k-1)*exp(-ix$theta^-1*q_1)) / (gamma(ix$k)^r_1 * ix$theta^{ix$k*s_1})
ix$y1<-ix$y.un1/sum(ix$y.un1, na.rm=T)

#The posterior probability at each grid location approximating p(k,theta|x1)
ix$y.un2<-(p_2^(ix$k-1)*exp(-ix$theta^-1*q_2)) / (gamma(ix$k)^r_2 * ix$theta^{ix$k*s_2})
ix$y2<-ix$y.un2/sum(ix$y.un2, na.rm=T)

#Plot the prior
ggplot(ix, aes(x=k, y=theta, fill=y.prior)) +
  geom_raster() +
  labs(x='Shape (k)',
       y='Scale (theta)',
       fill = 'Probability') +
  scale_fill_gradient(
    low = "#FFFFFF",
    high = "#0033FF",
    space = "Lab",
    na.value = "grey50",
    guide = "colourbar",
    aesthetics = "fill"
  )

# #Plot posterior 1
# ggplot(ix, aes(x=k, y=theta, fill=y1)) +
#   geom_raster() +
#   labs(x='Shape (k)',
#        y='Scale (theta)',
#        fill = 'Probability') +
#   scale_fill_gradient(
#     low = "#FFFFFF",
#     high = "#0033FF",
#     space = "Lab",
#     na.value = "grey50",
#     guide = "colourbar",
#     aesthetics = "fill"
#   )

# #Plot posterior 2
# ggplot(ix, aes(x=k, y=theta, fill=y2)) +
#   geom_raster() +
#   labs(x='Shape (k)',
#        y='Scale (theta)',
#        fill = 'Probability') +
#   scale_fill_gradient(
#     low = "#FFFFFF",
#     high = "#0033FF",
#     space = "Lab",
#     na.value = "grey50",
#     guide = "colourbar",
#     aesthetics = "fill"
#   )

# ix$beta=1/ix$theta
#Plot posterior 1
p1=ggplot(ix, aes(x=k, y=theta, fill=y1))+
  geom_tile() +
  labs(x=as.expression(bquote(alpha)),
       y=as.expression(bquote(1/beta)),
       fill = 'Posterior \nProbability') +
  scale_fill_gradient(
    limits=c(0, 0.012),
    breaks=c(0,0.01),
    low = "white",
    high = "black",
    space = "Lab",
    na.value = "grey50",
    guide = "colourbar",
    aesthetics = "fill"
  )+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        # legend.position = "none",
        axis.title = element_text(face="bold"),
        axis.text=element_text(size=12))+
  ggtitle("Five Observations of 20")

#Plot posterior 2
p2=ggplot(ix, aes(x=k, y=theta, fill=y2)) +
  geom_raster() +
  labs(x=as.expression(bquote(alpha)),
       y=as.expression(bquote(1/beta)),
       fill = 'Posterior \nProbability') +
  scale_fill_gradient(
    limits=c(0, 0.012),
    breaks=c(0,0.01),
    low = "white",
    high = "black",
    space = "Lab",
    na.value = "grey50",
    guide = "colourbar",
    aesthetics = "fill"
  )+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        # legend.position = "bottom",
        axis.title = element_text(face="bold"),
        axis.text=element_text(size=12))+
  ggtitle("Five Observations of 100")

p1+p2+plot_layout(ncol = 1)

ggsave(file="f_gammaheat.pdf",width = 4,height = 5)

#Now let's look a the 'posterior predictives'
#that is, the likelihood of new data x_prime, given the different posteriors on k and theta
df.pred<-data.frame(x = seq(0,250, length.out = 300))

#Make a big matrix to store all the individual parameter likelihood functions in across a plotting range of x_primes
ys<-matrix(NA, nrow = nrow(ix), ncol=length(df.pred$x))
for (i in 1:nrow(ix))
{
  ys[i,]<-dgamma(df.pred$x, shape = ix$k[i], scale = ix$theta[i])
}

#Weigh each by its posterior probability and sum them up to get the marginal
df.pred$five_obs_of_20<-colSums(sweep(ys, 1,  ix$y1, '*'))
#Separately for the case exposed to larger observations
df.pred$five_obs_of_100<-colSums(sweep(ys, 1,  ix$y2, '*'))

#Merge into one dataframe
df.pred.l<-df.pred %>% gather(key, y, five_obs_of_20:five_obs_of_100)

#Create an auxilliary dataframe showing some faint examples of realizations of the marginal
#by sampling gammas according to their posterior probability
par_ix1<-sample(1:nrow(ix), 100, p=ix$y1, replace = T)
par_ix2<-sample(1:nrow(ix), 100, p=ix$y2, replace = T)
df.aux<-data.frame(x = df.pred$x)
for (i in 1:100)
{
  df.aux[[paste0('y1_',i)]]<-dgamma(df.pred$x, shape = ix$k[par_ix1[i]], scale = ix$theta[par_ix1[i]])
}
#Do it again for the posterior based on the larger dataset
for (i in 1:100)
{
  df.aux[[paste0('y2_',i)]]<-dgamma(df.pred$x, shape = ix$k[par_ix2[i]], scale = ix$theta[par_ix2[i]])
}

#Merge into one very long auxilliary dataframe
df.aux.l<-df.aux %>% gather(y_key, five_obs_of_20, y1_1:y1_100) %>% select(-contains('y2'))
df.aux.tmp<-df.aux %>% gather(y_key, five_obs_of_100, y2_1:y2_100)
df.aux.l$five_obs_of_100<-df.aux.tmp$five_obs_of_100
df.aux.ll <- df.aux.l %>% gather(key, y, five_obs_of_20:five_obs_of_100)
df.aux.ll$y_key[df.aux.ll$key=='five_obs_of_100']<-rep(paste0('y2_',1:100), each = 300)
head(df.aux.ll)

#Plot everything together
ggplot(df.pred.l, aes(x=x, y=y, colour = key)) +
  geom_line(data = df.aux.ll, aes(y=y, group = y_key, colour = key), alpha = 0.05) +
  geom_line(size = 1, aes(group = key), colour = 'black') +
  geom_line(size = 1, linetype = 'dashed') +
  coord_cartesian(ylim = c(0,0.1)) +
  theme_bw()


df.pred.l=df.pred.l %>% 
  mutate(key=factor(key,levels = c("five_obs_of_20","five_obs_of_100"),
                    labels=c("Five Observations of 20","Five Observations of 100")))
df.aux.ll=df.aux.ll %>% 
  mutate(key=factor(key,levels = c("five_obs_of_20","five_obs_of_100"),
                    labels=c("Five Observations of 20","Five Observations of 100")))
#Plot everything together

ggplot(df.pred.l, aes(x=x, y=y, colour = key))+
  geom_line(data = df.aux.ll, aes(y=y, group = y_key, colour = key),linetype="solid", alpha = 0.03) +
  geom_line(size = 1) +
  geom_line(size = 1, aes(group = key), colour = 'black',linetype="dashed") +
  theme_classic()+
  theme(legend.position = c(0.7,0.8),
        axis.text=element_text(size=12),
        axis.title = element_text(face="bold"),
        legend.title=element_blank(),
        legend.background = element_rect(fill='transparent'),
        legend.text.align = 0)+
  xlab("Time")+
  ylab("Predicted Density")+
  scale_color_manual(values=alpha(c("#4DAF4A","#FF7F00"),0.9))
  # scale_linetype_manual(values=c("solid", "twodash", "dashed", "dotted"))
  # theme_bw()

ggsave(file="f_gammaland.pdf",width = 4.5,height = 3.5)
