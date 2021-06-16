import sage.all
from sage.rings.real_mpfr import RR
from sage.misc.functional import numerical_approx
from ore_algebra import *
import sys
sys.path.append("/home/agoyer/Work/MyPrograms")
from diffop_factorization import *
from diffop_factorization.examples import fcc3
from ore_algebra.examples import fcc
import cProfile
import pstats

fcc = [0, 0, 0, fcc3, fcc.dop4, fcc.dop5, fcc.dop6]

operators = [(fcc[3], 'fcc3'), \
(fcc[4], 'fcc4'), \
(fcc[5], 'fcc5'), \
(fcc[6], 'fcc6'), \
(fcc[3]**2, 'fcc3^2'), \
(fcc[4]*fcc[3], 'fcc4*fcc3'), \
(fcc[3]*fcc[4], 'fcc3*fcc4'), \
(fcc[4]**2, 'fcc4^2'), \
(fcc[3].lclm(fcc[4]), 'lclm(fcc3, fcc4)') ]

N = 5                       # number of loops
info = '21.06.16.16.53'     # Year.Month.Day.Hour.Min
o = open('sage.fcc.'+info+'.Nloops='+str(N), 'w')

key_tot = ('~', 0, '<built-in method builtins.exec>')
key_mono = None

fails = []

for (L, name) in operators[:2]:
    try:
        result = {'op': name}
        for cpt in range(N):
            cProfile.run('f = dfactor(L)', 'tmp_infos')
            s = pstats.Stats('tmp_infos')
            if key_mono is None:
                for key in s.stats.keys():
                    if key[2] == '_monodromy_matrices': key_mono = key
            tt , mt = s.stats[key_tot][3], s.stats[key_mono][3]
            if 'tot' not in result.keys() or tt < result['tot']:
                result['tot'] = numerical_approx(tt, digits=3)
                result['mono'] = numerical_approx(mt, digits=3)
                result['orders'] = [x.order() for x in f]
        o.write('operator: ' + result['op'] + '\n')
        o.write('total_time: ' + str(result['tot']) + '\n')
        o.write('mono_time: ' + str(result['mono']) + '\n')
        o.write('orders: ' + str(result['orders']) + '\n' + '\n')
    except:
        fails.append(name)

if not fails == []: print('fails:'+str(fails))

o.close()
quit
