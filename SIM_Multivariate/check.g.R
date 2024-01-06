

#Just to check the symmetry of g issue.
#Issue: we all agree g should be symmetric. However, Omega and mu are not, and g = Omega %*% mu %*% Omega'

#Rule: don't transpose when going against the direction of the arrow. Do when going direction of arrow



#Assume that Omega is full (I'm still not sure if this is right), but that mu is symmetric (as I've argued) and so is g.

#Then g = Omega' mu Omega, so let's see

(g <- matrix(c(.3, .1,
              .1, .2),nrow=2,byrow=T))

(Omega <- matrix(c(.15, .05,
                   .01, .10),nrow=2,byrow=T))

(mu <- matrix(c(.3, .1,
                .1, .3),nrow=2,byrow=T))

t(Omega) %*% mu %*% Omega



#Other params

(g <- matrix(c(.3, .1,
               .1, .2),nrow=2,byrow=T))

(Omega <- matrix(c(.5, .25,
                   .01, .10),nrow=2,byrow=T))

(mu <- matrix(c(.3, .1,
                .1, .3),nrow=2,byrow=T))

t(Omega) %*% mu %*% Omega