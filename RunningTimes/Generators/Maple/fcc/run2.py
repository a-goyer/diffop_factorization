N = 5                         # number of loops
info = '21.05.13.H.M'       # Year.Month.Day.Hour.Min

o = open('maple.fcc.'+info+'.Nloops='+str(N), 'w')

for i in range(3,6):
    t = Infinity
    for j in range(N):
        f = open('tmp/o'+str(j), 'r')
        l = f.readlines()
        for u in l:
            if u[:8] == '### fcc'+str(i): t = min(t, RR(u[13:-1]))
    t = numerical_approx(t, digits=3)
    o.write('fcc'+str(i)+':\n')
    o.write(str(t)+'\n')

o.close()
quit
