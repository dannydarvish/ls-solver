# for E in range(300,1501,10):
# for E in range(300,1010,10):
# for E in range(900,1000,1):
for E in range(750,1201,10):
    if E == 700:
        E = 690

    lines = '''#!/bin/bash
#PBS -N Lippmann-Schwinger-inf-vol
#PBS -q red
#PBS -l walltime=48:00:00
#PBS -l cput=256:00:00
#PBS -l nodes=1:ppn=8
#PBS -j oe
#PBS -o /latticeQCD/raid6/darvish/ls-solver/work/jobs/E_%d.log
#PBS -e /latticeQCD/raid6/darvish/ls-solver/work/jobs/E_%d_err.log


''' % (E,E)
    lines += '''export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib:/home/darvish/.local/lib:/home/darvish/.local/lib64:/home/darvish/.local/libexec
export PATH=$PATH:/usr/local/bin

'''
    lines += 'cd /latticeQCD/raid6/darvish/ls-solver/work && /latticeQCD/raid6/darvish/ls-solver/work/test-1b2c %f' % E
    lines += '\n'

    f = open('E_%d_job' % E, 'w')
    f.write(lines)
    f.close()