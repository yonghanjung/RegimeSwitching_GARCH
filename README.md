# RegimeSwitching_GARCH

It is written in MATLAB. 
The original matlabcode can be found in 
"http://www.thomaschuffart.fr/?page_id=12"

Howver, as there are many errors in the code uploaded at the above link,
I modified the code for the problem of 
* Two regime 
* GARCH(1,1) 

<hr> 
For estimating the parameters of Regime switching GARCH(1,1), 
you have to run 
<code> [thetahat results struct]= swgarchest(data,flag,ORDERS,reg) </code>

* flag = 2 is recommended, 
* I only tested it for the 'reg = 2'
* I only tested it for ORDERS = [1,1,1]
