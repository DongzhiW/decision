job.java

Change the state from Employed to Unemployed or retired
```
m1 = ParameterInputs.EmpToUnempConst
					  + ParameterInputs.EmpToUnempAgeCoeff * Age
					  + ParameterInputs.EmpToUnempRateCoeff * unemploymentRate
					  + ParameterInputs.EmpToUnempIncomeCoeff * this.income
					  + ParameterInputs.EmpToUnempHealthCoeff * (1 - HealthFactor)
					  + ParameterInputs.EmpToUnempEduCoeff * educ;
			  m1 = 1.0 - 1.0 / (1.0 + Math.exp(m1));

			  m2 = ParameterInputs.EmptoRetconst +
						  ParameterInputs.EmptoRetAgeCoeff * Age +
						  ParameterInputs.EmptoRetHHCoeff * hhSize +
						  ParameterInputs.EmptoRetEduCoeff * educ +
						  ParameterInputs.EmptoRetHealthCoeff * (1 - HealthFactor) +
						  ParameterInputs.EmptoRetAssetsCoeff * (assets - debt) +
						  ParameterInputs.EmptoRetExpenseCoeff * expense;
			  m2 = 1.0 - 1.0 / (1.0 + Math.exp(m2));
```

VA.java

```
int Cost = (int)((ParameterInputs.VALapseprob * lapsefv + (1.0 - ParameterInputs.VALapseprob)*DBfv)* ParameterInputs.VAannCushion);

if(d1==0) rndm=ParameterInputs.VAattention_t1;
	else if(d1==1)rndm=ParameterInputs.VAattention_t2;
	else if(d1==2)rndm=ParameterInputs.VAattention_t3;
	rndm *= attentionMult;


probAccept = ParameterInputs.VABuyOutcoeff1
			 + ParameterInputs.VABuyOutcoeff1 * dbTerm
			 + ParameterInputs.VABuyOutcoeff3 * savingsPerformance
			 + ParameterInputs.VABuyOutcoeff4 * ibTime
			 + ParameterInputs.VABuyOutcoeff5 * offerRt
			 + ParameterInputs.VABuyOutcoeff6 * vsPrior;

	 probAccept = 1.0 / (1.0 + Math.exp(probAccept));

	 if(probAccept < draw) return true;
	 else return false;	 

```

```
public boolean regularWDAttention(int age, Date dt1, double durUnemp){
 	 double chance = rand.nextDouble();
	 double prbAtten = ParameterInputs.regWDAttenIntercept
			+ ParameterInputs.regWDAttenAge * (age * 1.0)
			+ Math.pow(ParameterInputs.regWDAttenAnn,  Commonmethods.getMonthsSinceLastAnniversary(EffDate, dt1))
	 		+ ParameterInputs.regWDAttenUnempl * durUnemp;
 	prbAtten = 1.0 / (1.0 + Math.exp(-1.0 * prbAtten));
 	if(chance < prbAtten){ return true; }
 	else {return false;}
 }
 ```
 
 ```
 double probRegWD = ParameterInputs.regWdDec_intercept
        + ParameterInputs.regWdDec_time * isCorrectTime(mf, age)
        + ParameterInputs.regWdDec_db * dbFactor
        + ParameterInputs.regWdDec_tax * taxStatus
        + ParameterInputs.regWdDec_income * incomeStrain;

     probRegWD = 1.0 / (1.0 + Math.exp(-1.0*probRegWD));
     if(rand.nextDouble() < probRegWD) { return true;}
     else { return false;}
```

Person.java

```
double probD4D = ParameterInputs.wdh_t1Intercept
				+ ParameterInputs.VAwdHxWeight * v.ifWDTaken()
				+ ParameterInputs.VAmktPerfWeight * (v.getStanding() - 1.0)
				+ ParameterInputs.VACashNeedageCoef * this.age
				+ ParameterInputs.CreditBurnedCoef * Commonmethods.zidz(d.total_outstanding_debt(), d.maximum_loan_capacity());
		probD4D = 1.0/(1.0+Math.exp(-1.0*probD4D));
```

```
private boolean t2FundingDecision(double args[]){
		count[12][0]++;
		double x = ParameterInputs.wdh_t2Intercept
				+ args[0] * ParameterInputs.RetWhAgeCoef 			// 401k: age dependency
				+ args[1] * ParameterInputs.RetWhTotNCoef 			// 401k: sqrt((shortfall / 401k max wd) - 1)
				+ args[2] * ParameterInputs.RetWithdrawalEfficiency // 401k: realized withdrawal / withdrawn amt
				+ args[3] * ParameterInputs.VAEffWhCoef 			// VA: realized withdrawal / withdrawn amt
				+ args[4] * ParameterInputs.VAnlgWhCoef 			// VA: won't break NLG or LAPSE
				+ args[5] * ParameterInputs.VAneedPctCoef 			// VA: sqrt((shortfall / max pro-rata wd) - 1)
				+ args[6] * ParameterInputs.VAdbCoef 				// VA: DB importance
				+ args[7] * ParameterInputs.VAitmCoef;				// VA: moneyness

		 x = 1.0/(1.0+Math.exp(-1.0*x));
 		 if(rand.nextDouble() < x) 
			 return true;  // 401k Withdrawal
		 else { 
			 count[12][1]++;
			 return false; // VA Pro-Rata Withdrawal
		 } 
	}
```

```
double x = ParameterInputs.vaSurrAttenConst + ParameterInputs.VASurrChancoeff * va.getChannelObject().isWholesale();
	 x = 1.0/(1.0+Math.exp(-1.0*x));
	 x *= Math.pow(ParameterInputs.vaSurrAttenAnn, monthsAfterAnniversary);
```     

voluntarily lapse/given the market growth rate is much higher than the rollup rate
```
double xlapse = ParameterInputs.SurrAttentionConst
					+ (av/100000)*ParameterInputs.AVAttentionCoeff
		            + wltshr*ParameterInputs.AssetAttentionCoeff
		            + ParameterInputs.VASurrRUcoeff*(va.BB[1].getRollupRate() - m.y.gettenYearRate())
		            + ParameterInputs.VASurrRatiocoeff * itm
		            + ParameterInputs.VASurrDBcoeff * b
		            + ParameterInputs.VASurrCWCcoeff * ShareClassFees.getCWC(va.getSCName(), va.getDuration(d))
		            + ParameterInputs.VASurrCWCDiffcoeff * CWCdiff;
		    
		    return  (1.0/(1.0+Math.exp(-1.0*xlapse)) > a);
 ```
