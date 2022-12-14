from scipy import optimize

class Make_Up_Flow_Balaner(object):
    def __init__(self) -> None:
        pass

    def solve(self,target_CaO_conversion):
        # target_CaO_conversion should be less than 0.66, otherwise no result
        X=target_CaO_conversion
        r=optimize.brentq(self._f,0.0001,1,args=(X,),rtol=1e-4)

        results={}
        results["make_up_ratio"]=r
        # results["calcination_fraction"]=X/(r+(1-r)*X) #fcalc
        results["calcination_fraction"]=1
        # results["caco3_mole_fraction_calciner_outlet"]=r*(1-X) #Y
        results["caco3_mole_fraction_calciner_outlet"]=0 #Y
        return results


    def _f(self,make_up_ratio,target_CaO_conversion):
        # correlation between average CaO conversion and make_up ratio
        # reference: Chemical Engineering Journal 156(2010):388-394
        a1=0.1045
        f1=0.9822
        a2=0.7786
        f2=0.7905
        b=0.07709
        r=make_up_ratio
        X=0.1
        X_1=0
        eps=1e-3
        while abs(X-X_1)>eps:
            X_1=X
            # fcalc=X/(r+(1-r)*X)
            fcalc=1
            r0=r*(1-fcalc)/(r+(1-r)*fcalc)
            X=(r+(1-r)*r0)*fcalc*(a1*f1*f1/(r+(1-r)*fcalc*(1-f1))+a2*f2*f2/(r+(1-r)*fcalc*(1-f2))+b/r)
        return target_CaO_conversion-X
