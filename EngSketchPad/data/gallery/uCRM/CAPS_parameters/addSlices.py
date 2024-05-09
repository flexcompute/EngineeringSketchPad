import numpy as np

baselineDesign = {'kulfan8x8_inboard102_aupper4':0.16742,
                  'kulfan8x8_inboard102_adelta8':-0.09452,
                  'kulfan8x8_inboard102_aupper1':0.11876,
                  'kulfan8x8_inboard102_adelta6':-0.28262,
                  'kulfan8x8_inboard102_aupper2':0.09047,
                  'kulfan8x8_inboard102_adelta4':-0.01000,
                  'kulfan8x8_inboard102_aupper7':0.10596,
                  'kulfan8x8_inboard102_aupper5':0.07238,
                  'kulfan8x8_inboard102_adelta7':-0.16159,
                  'kulfan8x8_inboard102_aupper3':0.12322,
                  'kulfan8x8_inboard102_adelta5':-0.48747,
                  'kulfan8x8_inboard102_aupper8':0.34921,
                  'kulfan8x8_inboard102_adelta3':-0.37575,
                  'kulfan8x8_inboard102_aupper6':0.23458,
                  'kulfan8x8_inboard102_adelta1':-0.22402,
                  'kulfan8x8_inboard102_adelta2':-0.22501,
                  'kulfan8x8_outboard152_aupper4':0.06044,
                  'kulfan8x8_outboard152_adelta8':-0.15837,
                  'kulfan8x8_outboard152_aupper1':0.08339,
                  'kulfan8x8_outboard152_adelta6':-0.23355,
                  'kulfan8x8_outboard152_aupper2':0.05207,
                  'kulfan8x8_outboard152_adelta4':-0.01000,
                  'kulfan8x8_outboard152_aupper7':0.06046,
                  'kulfan8x8_outboard152_aupper5':0.05000,
                  'kulfan8x8_outboard152_adelta7':-0.15823,
                  'kulfan8x8_outboard152_aupper3':0.08085,
                  'kulfan8x8_outboard152_adelta5':-0.29372,
                  'kulfan8x8_outboard152_aupper8':0.10464,
                  'kulfan8x8_outboard152_adelta3':-0.22624,
                  'kulfan8x8_outboard152_aupper6':0.16611,
                  'kulfan8x8_outboard152_adelta1':-0.19149,
                  'kulfan8x8_outboard152_adelta2':-0.07113,
                  'kulfan8x8_inboard102_ztail2':0.00500,
                  'kulfan8x8_outboard152_ztail2':-0.00135}

Nsect = 12 # number of slices, including beginning and end

name1 = ['kulfan8x8_inboard102','kulfan8x8_outboard152']
name2 = ['aupper1','aupper2','aupper3','aupper4','aupper5','aupper6','aupper7','aupper8',
         'adelta1','adelta2','adelta3','adelta4','adelta5','adelta6','adelta7','adelta8','ztail2']

Y = np.linspace(0.,1.,Nsect)
Ystr = ''
TWISTstr = ''
TWISTdespmtr = ''
for i in range(len(Y)-2):
    Ystr += ('%f; ' % Y[i+1])
    TWISTstr += ('Wing:outboardSec:twistInterp%.0f; ' % (i+1))
    TWISTdespmtr += ("""
DESPMTR Wing:outboardSec:twistInterp%.0f %f""" % (i+1,0.837972+Y[i+1]*(4.612053-0.837972)))
    for j in range(len(name2)):
        newDVname = ('kulfan8x8_interp%.0f_%s' % (i+1,name2[j]))
        newDVval = baselineDesign[('%s_%s' % (name1[0],name2[j]))] + \
                   Y[i+1]*(baselineDesign[('%s_%s' % (name1[1],name2[j]))]-baselineDesign[('%s_%s' % (name1[0],name2[j]))])
        print('despmtr %s %f' % (newDVname,newDVval))
print('\n\n')
print(Ystr)
print(TWISTstr)
print(TWISTdespmtr)
print('\n\n')

for i in range(len(Y)-2):
    for j in range(len(name2)):
        newDVname = ('kulfan8x8_interp%.0f_%s' % (i+1,name2[j]))
        index = name2[j][-1]
        print('set Wing:outboardSec:%s[%s,%.0f] %s' % (name2[j][:-1],index,i+4,newDVname))

