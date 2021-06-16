from ore_algebra import *
from diffop_factorization import *
from ore_algebra.examples import ssw

key_tot = ('~', 0, '<built-in method builtins.exec>')
key_mono = ('monodromy.py', 406, 'monodromy_matrices')

N = 5                         # number of loops
info = '21.05.13.H.M'       # Year.Month.Day.Hour.Min
o = open('sage.ssw.'+info+'.Nloops='+str(N), 'w')

for i in range(2):
    for j in range(2):
        for k in range(1,20):
            print('k i j =', k, i, j)
            try:
                L = ssw.dop[k,i,j]
                result = {'op': 'ssw.dop[' + str(k) + ',' + str(i) + ',' + str(j) + ']'}
                mono = False
                for cpt in range(N):
                    s = %prun -rq f = dfactor(L)
                    if key_mono in s.stats.keys() :
                        tt , mt = s.stats[key_tot][3], s.stats[key_mono][3]
                        if 'tot' not in result.keys() or tt < result['tot']:
                            result['tot'] = numerical_approx(tt, digits=3)
                            result['mono'] = numerical_approx(mt, digits=3)
                            result['orders'] = [x.order() for x in f]
                            mono = True
                    else:
                        tt = s.stats[key_tot][3]
                        if 'tot' not in result.keys() or tt < result['tot']:
                            result['tot'] = numerical_approx(tt, digits=3)
                            result['orders'] = [x.order() for x in f]
                            mono = False
                o.write('operator: ' + result['op'] + '\n')
                o.write('total_time: ' + str(result['tot']) + '\n')
                if mono:
                    o.write('mono_time: ' + str(result['mono']) + '\n')
                else:
                    o.write('mono_time: None' + '\n')
                o.write('orders: ' + str(result['orders']) + '\n' + '\n')
            except Exception:
                print('!!!!! ssw.dop['+str(k)+','+str(i)+','+str(j)+'] failed. !!!!!')

o.close()
quit
