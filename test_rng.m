rng('default')
s = rng; % save the current random generator settings
x = rand(1,5)
xnew = rand(1,5)

rng(s) % restore the original random generator settings based on s
xold = rand(1,5)

t = rng('default'); % save the current random generator settings, then set the random generator settings to 'default'
y = rand(1,5)
ynew = rand(1,5)

rng(t) % restore the random generator settings based on t
xnext = rand(1,5)