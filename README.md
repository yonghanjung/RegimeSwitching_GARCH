# RegimeSwitching_GARCH

It is written in MATLAB. 
The original matlab code that I refered can be found in "http://www.thomaschuffart.fr/?page_id=12",
and it is the same code used by Juri Marcucci, for "Forecasting Stock Market Volatility with Regime-Switching GARCH Models". (M. Juri. Forecasting stock market volatility with regime-switching garch models. Studies in Nonlinear Dynamics & Econometrics, 9(4), 2005) 

As there are many errors in the code (Not runnable) I debugged and modified the code for the case of 
* two regimes 
* GARCH(1,1) 

<hr> 
For estimating the parameters of Regime switching GARCH(1,1), 
you have to run 
<code> [thetahat results struct]= swgarchest(data,flag,ORDERS,reg) </code>

* flag = 2 is recommended, 
* I only tested it for the 'reg = 2'
* I only tested it for ORDERS = [1,1,1] or  [0,1,1]
  * [1,1,1] means GARCH(1,1) with the constant (mu) existed, [0,1,1] means GARCH(1,1) with the constant (mu) NOT existed

<hr> 
As the code is SUPER slow, I will port it into the Python code. 
