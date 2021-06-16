N = 5                         # number of loops
info = '21.04.26.H.M'       # Year.Month.Day.Hour.Min

o = open('maple.ssw.'+info+'.Nloops='+str(N), 'w')

for i in range(2):
    for j in range(2):
        for k in range(1,20):
            t = 0
            for cpt in range(N):
                f = open('tmp/o'+str(cpt), 'r')
                l = f.readlines()
                for u in l:
                    s = '##### dop['+str(k)+','+str(i)+','+str(j)+']'
                    if u[:16] == s[:16]: t += RR(u[23:-1])
            t = (t/N).n(20)
            o.write('dop['+str(k)+','+str(i)+','+str(j)+']:\n')
            o.write(str(t)+'\n')

o.close()
quit
