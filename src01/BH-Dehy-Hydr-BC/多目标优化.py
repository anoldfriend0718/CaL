import numpy as np
import geatpy as ea
import matplotlib.pyplot as plt
plt.rcParams['font.sans-serif'] = ['DejaVu Sans']

class MyProblem(ea.Problem):
    def __init__(self):
        name = 'BNH'
        M=2
        maxormins = [1]*M
        Dim = 2
        varTypes = [0]*Dim
        lb = [0]*Dim
        ub = [5,3]
        lbin = [1]*Dim
        ubin = [1]*Dim
        ea.Problem.__init__(self,name,M,maxormins,Dim,varTypes,lb,ub,lbin,ubin)
    
    def aimFunc(self,pop):
        Vars =pop.Phen
        x1 = Vars[:,[0]]
        x2 = Vars[:,[1]]
        f1 = 4*x1**2+4*x2**2
        f2 = (x1-5)**2+(x2-5)**2
        pop.CV = np.hstack([(x1-5)**2+x2**2-25,-(x1-8)**2-(x2-3)**2+7.7])
        pop.ObjV = np.hstack([f1,f2])

    def calReferObjV(self):
        N = 10000
        x1 = np.linspace(0,5,N)
        x2 = x1.copy()
        x2[x1>=3]= 3
        return np.vstack((4*x1**2+4*x2**2,(x1-5)**2+(x2-5)**2)).T
    
problem = MyProblem()
Encoding = 'RI'
NIND = 100
Field = ea.crtfld(Encoding,problem.varTypes,problem.ranges,problem.borders)
population = ea.Population(Encoding,Field,NIND)
myAlgorithm = ea.moea_NSGA2_templet(problem,population)
myAlgorithm.mutOper.Pm = 0.2
myAlgorithm.recOper.XOVR = 0.9
myAlgorithm.MAXGEN = 200
myAlgorithm.logTras = 1
myAlgorithm.verbose = False
myAlgorithm.drawing = 1
[NDSet,population] = myAlgorithm.run()
NDSet.save()
print('用时 %s 秒' % myAlgorithm.passTime)
#print('非支配个体数：%d 个' % NDSet.sizes) if NDSet.sizes! = 0 else print('没有找到可行解！')
if myAlgorithm.log is not None and NDSet.sizes !=0:
    print('GD',myAlgorithm.log['gd'][-1])
    print('IGD',myAlgorithm.log['igd'][-1])
    print('HV',myAlgorithm.log['hv'][-1])
    print('Spacing',myAlgorithm.log['spacing'][-1])
metricName = [['igd'],['hv']]
Metrics = np.array([myAlgorithm.log[metricName[i][0]] for i in range(len(metricName))]).T
ea.trcplot(Metrics,labels=metricName,titles = metricName)

#import matplotlib.font_manager as fm
 
# 获取所有可用的字体
#fonts = fm.findSystemFonts(fontpaths=None, fontext='ttf')
 
# 打印所有字体的名称和路径
#for font in fonts:
#    prop = fm.FontProperties(fname=font)
#    print(f"Font Name: {prop.get_name()}, Font Path: {font}")