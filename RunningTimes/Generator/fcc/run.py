from ore_algebra import *
from diffop_factorization import *
from diffop_factorization.examples import fcc3
from ore_algebra.examples import fcc

fcc = [0, 0, 0, fcc3, fcc.dop4, fcc.dop5, fcc.dop6, fcc3*fcc3, fcc.dop4*fcc3, fcc3*fcc.dop4, fcc.dop4*fcc.dop4]

N = 1                       # number of loops
info = '21.05.13.H.M'     # Year.Month.Day.Hour.Min
o = open('sage.fcc.'+info+'.Nloops='+str(N), 'w')

key_tot = ('~', 0, '<built-in method builtins.exec>')
key_mono = ('monodromy.py', 406, 'monodromy_matrices')

fails = {}

for i in range(3, 11):
    print('i =', i)
    L = fcc[i]
    if i<7:
        result = {'op': 'fcc'+str(i)}
    elif i==7:
        result = {'op': 'fcc3**2'}
    elif i==8:
        result = {'op': 'fcc4*fcc3'}
    elif i==9:
        result = {'op': 'fcc3*fcc4'}
    elif i==10:
        result = {'op': 'fcc4**2'}
    for cpt in range(N):
        s = %prun -rq f = dfactor(L)
        tt , mt = s.stats[key_tot][3], s.stats[key_mono][3]
        if 'tot' not in result.keys() or tt < result['tot']:
            result['tot'] = numerical_approx(tt, digits=3)
            result['mono'] = numerical_approx(mt, digits=3)
            result['orders'] = [x.order() for x in f]
    o.write('operator: ' + result['op'] + '\n')
    o.write('total_time: ' + str(result['tot']) + '\n')
    o.write('mono_time: ' + str(result['mono']) + '\n')
    o.write('orders: ' + str(result['orders']) + '\n' + '\n')

o.close()
quit
