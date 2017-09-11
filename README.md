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
 public boolean annuitizationDecision(int Age, ActuarialFunc AF, boolean gender){
	 boolean annuitize = false;
	 int perceptionBias = rand.nextInt(6); //how long does the policyholder think they will live for sure
	 if(AF.getHealthStatus() == HealthState.Poor) perceptionBias = 0;
	 AF.setDiscAdj_t(ParameterInputs.VAHyperDiscFactor); //override discount factor in AF so we have hyperbolic discounting
	 int Value = (int) AF.calc_nCertLifeAnn(Age, perceptionBias) * this.calcAnnPayment(AF, gender, Age); //expected pv of annuity times payment amount
	 int DBfv = 0;

	 if(BB[BenefitType.DB.getIndex()].RiderType.equals("MAX"))
	 {//RU BB grows so we need to split into an increasing term life and a deferred term life with a higher benefit
		 int dbGrowPeriod = Math.max(BB[BenefitType.IB.getIndex()].MaxAge - Age, 0); //DB Rolls up to end of IB
		 double dbInc =  BB[BenefitType.DB.getIndex()].BB[0] * AF.calcIncLife(Age, dbGrowPeriod, BB[BenefitType.DB.getIndex()].RollupRate); //increasing life ins from growing db
		 double dbTerm = (BB[BenefitType.DB.getIndex()].BB[0] * Math.pow(1+BB[BenefitType.DB.getIndex()].RollupRate, dbGrowPeriod)); // the DB has rolled up to the max ib age
         dbTerm *= AF.calc_defTermLife(Age, dbGrowPeriod, Math.max(BB[BenefitType.DB.getIndex()].MaxAge - BB[BenefitType.IB.getIndex()].MaxAge, 0)); // term life ending at max db age with ru benefit
         DBfv = (int) (dbInc + dbTerm);
	 }
	 else
     {//ROP and Ratchet dont grow so we can treat the whole benefit as one longer term policy
         int t = Math.max(BB[BenefitType.DB.getIndex()].MaxAge - Age, 0);
         DBfv = (int) (Math.max(BB[BenefitType.DB.getIndex()].BB[1],BB[BenefitType.DB.getIndex()].BB[2]) * AF.calcLife(Age, t));
     }

	 int lapsefv = Commonmethods.sumarray(AccountValue) ;
	 int Cost = (int)((ParameterInputs.VALapseprob * lapsefv + (1.0 - ParameterInputs.VALapseprob)*DBfv)* ParameterInputs.VAannCushion);
	 //@TODO Determine how comparable value and cost actually are
	 if(Value > Cost) annuitize = true;
	 AF.setDiscAdj_t(1); //reset the hyperbolic discount to 1 - rational tvm
	 return annuitize;
 }

 public boolean buyoutDecision(int Age, HealthState h, Date d, boolean gender, int total_assets, double netcashflow, int PersonId, MarketFactors m, ActuarialFunc AF){
	 double offerAmt = updateOfferAmount(PersonId, m, d); 
	 this.buyoutelig = false;
	 this.buyoutamt = 0;
	 boolean AcceptBuyout = false;
	 double anchorRef = 1.0; //not used unless has offer anchor is true
     if(this.hasOfferAnchor){ anchorRef = this.offerAnchor;}
	 if (offerAmt > 0){
		 double rndm = 0;
		 this.buyoutelig = true;
		 this.buyoutamt = (int) offerAmt;
		 double attentionMult = 1.0;
		 if(hasOfferAnchor) //the policyholder will pay less attention to subsequent offers, but attention recovers overtime
		 {
		 	attentionMult = (1.0 - Math.pow(ParameterInputs.VAattention_serial, this.yrsSinceLastOffer)); //reduction to attention for serial offers
		 }
		 int d1=Commonmethods.getDateDifferenceinCalMonths(OfferAmountLookup.getOfferDate(PersonId), d);
		 if(d1==0) rndm=ParameterInputs.VAattention_t1;
		 else if(d1==1)rndm=ParameterInputs.VAattention_t2;
		 else if(d1==2)rndm=ParameterInputs.VAattention_t3;
		 rndm *= attentionMult;
		 double a = rand.nextDouble();
		 if(a < rndm){
		 	AcceptBuyout = offerDecision(h, netcashflow, total_assets,Math.max(BB[BenefitType.IB.getIndex()].MaxAge - Age,0), anchorRef, offerAmt);
            newOfferAnchor =  offerAmt / (BB[BenefitType.IB.getIndex()].IBUseBB() * 1.0);
            updateOffer = true;
         }
         if((d1>=2)&&updateOffer){
		    offerAnchor = newOfferAnchor;
		    yrsSinceLastOffer = 0;
		    hasOfferAnchor = true;
            newOfferAnchor = 0;
            updateOffer = false;
         }
	 }
	 return AcceptBuyout;
 }

	/**
	 * Function to calculate the probability that the policyholder accepts the offer and simulation of decision
	 * <p>
	 *     The function is based on the original logic developed in conjunction with PWC. The policyholder considers
	 *     how they are feeling if they have a max db, how much time the have until the ib goes away, how their cash flow
	 *     has been over recent time, and the relative value of the offer compared to what they think the contract is worth
	 * </p>
	 * <p>
	 *     The probability takes the form 1/1+e(x) and the draw uses the random uniform distribution
	 * </p>
	 * @param h current health state
	 * @param netcashflow ewma of the momthly net cash flow for the ph
	 * @param assets total asset balances
	 * @param ibTime time left until IB expires
	 * @param anchorRatio offer / ib bb in last offer the ph saw
     * @param offerAmt the current buyout offer amount
	 * @return boolean indicating whether or not the policyholder has accepted the offer
	 * @todo confirm decision to have no offer ratio comparison term for first time offers
	 */
 private boolean offerDecision(HealthState h, double netcashflow, int assets, int ibTime, double anchorRatio, double offerAmt)
 {
	 int Health = 0;
	 int MAXDB = 0;
	 double probAccept;
	 double draw = rand.nextDouble();

	 if(h==HealthState.Poor) Health = 1;
	 if(BB[BenefitType.DB.getIndex()].RiderType.equals("MAX")) MAXDB = 1;

	 double dbTerm = MAXDB * Health; //dummy variable = 1 if sick and max db
	 double savingsPerformance = netcashflow/(1.0 * assets); //ratio of monthly net cash flow to total assets
	 double offerRt = offerAmt/(BB[BenefitType.IB.getIndex()].IBUseBB() * 1.0); //ratio of current offer amount to the IB benefit base
	 double vsPrior = 0.0; //holds the difference in the ratio between offer and IB BB in the current offer versus the last offer the ph noticed
	 if(this.hasOfferAnchor) vsPrior = offerRt - anchorRatio;

	 probAccept = ParameterInputs.VABuyOutcoeff1
			 + ParameterInputs.VABuyOutcoeff1 * dbTerm
			 + ParameterInputs.VABuyOutcoeff3 * savingsPerformance
			 + ParameterInputs.VABuyOutcoeff4 * ibTime
			 + ParameterInputs.VABuyOutcoeff5 * offerRt
			 + ParameterInputs.VABuyOutcoeff6 * vsPrior;

	 probAccept = 1.0 / (1.0 + Math.exp(probAccept));

	 if(probAccept < draw) return true;
	 else return false;
 }
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
