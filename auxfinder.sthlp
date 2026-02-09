{smcl}
{* 16jan2026}{...}
{hi:help auxfinder}
{hline}

{title:Title}

{pstd}{hi:auxfinder} {hline 2} Finding auxiliary variables for imputation models

{title:Syntax}

{p 8 15 2}
{cmd:auxfinder} {help varlist:{it:varlist}} {ifin}, {cmd:testvars}({help varlist:{it:varlist}}) [{cmd:}
{help auxfinder##comopt:{it:options}}]

{synoptset 25 tabbed}{...}
{marker comopt}{synopthdr:options}
{synoptline}
{synopt :{opt test:vars(varlist)}}List of variables tested as potential auxiliaries (required)
  {p_end}
{synopt :{opt corrm:iss(#)}}Required correlation between the test variable and the missing indicator of the target; default is {cmd: corrmiss(0.15)}
  {p_end}
{synopt :{opt corrv:alue(#)}}Required correlation between the test variable and the target; default is {cmd: corrvalue(0.30)}
  {p_end}
{synopt :{opt gain(#)}}Share of cases that gain extra information by variable; default is {cmd: gain(0.15)} 
  {p_end}
{synopt :{opt max:miss(#)}}Share of missingness in the testing variable allowed in relation to the target variable; default is {cmd: maxmiss(0.15)} 
  {p_end}
{synopt :{opt nonl:inear(#)}}Testing for nonlinear (quadratic) influences of continuous test variables
  {p_end}
{synopt :{opt cat:limit(#)}}Number of levels required to be regarded as a continuous test variable; default is {cmd: catlevel(7)} 
  {p_end}
{synopt :{opt sav:ing(filename, ...)}}save results to a file
  {p_end}
{synopt :{opt folds(#)}}number of folds for {help lasso:{it:lasso}} / {help elasticnet:{it:elasticnet}}; default is {cmd: folds(20)} 
  {p_end}
{synopt :{opt seed(#)}}seed for RNG; relevant for {help lasso:{it:lasso}} / {help elasticnet:{it:elasticnet}}
  {p_end}
{synopt :{opt skip:lasso}}{help lasso:{it:lasso}} / {help elasticnet:{it:elasticnet}} not computed
  {p_end}
{synopt :{opt elastic:net}}use {help elasticnet:{it:elasticnet}} instead of {help lasso:{it:lasso}}
  {p_end}
{synopt :{opt lassoopts(string)}}options passed to {help lasso:{it:lasso}} / {help elasticnet:{it:elasticnet}}
  {p_end}
{synopt :{opt nois:ily}}Show more output
  {p_end}
  
{title:Description}

{pstd}{cmd:auxfinder} selects relevant auxiliary variables
for imputation models. An auxiliary variable is a variable that is not included
in the analytical estimation model, but is added to the imputation model
to improve the generation of imputed values. Ideally, an auxiliary variable has no
(or only very little) missing values itself and is strongly correlated with the
variables that need to be imputed. {cmd:auxfinder} implements a two-step selection process:
first, missing shares and correlations between target variables (variables that need to be imputed)
are tested. If these results are satisfactory, a {help lasso:{it:lasso}} model is estimated.
Only variables will be recommended as auxiliary variables that pass both checks. The
command is built to work with large datasets, where selection of variables based on
theoretical assumptions becomes difficult.
A working paper describing {cmd:auxfinder} in detail is available; see Bittmann (2026).

{pstd}Just as {cmd:mi impute}, {cmd:auxfinder} only considers the system missing value (.)
as missings in the target variables. Extended missing values (.a .b .c) are
ignored and need to be manually recoded before using {cmd:auxfinder}.


{title:Options}

{phang}{opt testvars(varlist)} specifies the variables that are potential candidates as auxiliaries.
Factor variable notation is allowed (but omit the c. prefix; do not specifiy interactions).
To test all variables in the dataset, use the specifier _all.
To assess whether a variable is continuous or categorical, the number of levels is used, see the 
option catlimit(). At least one variable must be specified.

{phang}{opt corrmiss(#)} specifies the degree of correlation between the binary missingness indicator of the target
and the testing variable. A high correlation means that a testing variable is able to predict whether a case
has a missing value on the target variable or not, which benefits the imputation procedure.
Test variables with higher levels of correlation are preferred. To compute the correlation values, the square
roots of R² values are utilized.

{phang}{opt corrvalue(#)} specifies the degree of correlation between the target and the testing variable
to be eligible for selection. Variables with higher levels of correlation are preferred. To compute the
correlation values, the square roots of (pseudo) R² values are utilized. For continuous target variables,
linear (OLS) models are tested using {cmd:regress}. For categorical targets, {cmd:mlogit} is utilized.
While lowering the correlation will produce more potential candidates, their benefit can be lower.

{phang}{opt gain(#)} specifies the amount of additional information a testing variable needs to add for a given
target variable. The gain is the share of cases with missing values on the target variable that have a valid
value on the testing variable. Ideally, the gain is 1 (100%) and values closer to 1 are better.

{phang}{opt maxmiss(#)} Optimally, auxiliaries have no missing values. However, some amount of missingness in relation
to the target variable can be allowed. If this share is too high, the following lasso models might not converge due to a
limited number of information available. The default value of 0.15 means that a testing variable can have a share of missing
values that is 15% larger than the share of missing values in the target variable.

{phang}{opt nonlinear(#)} automatically tests for nonlinear influences of testing variables. While a correlation usually assumes a linear
association between two variables, quadratic influences are often present in real data, meaning that the associations
can be (inverted) U-shaped. This can be tested using regression models by introducing the squared term of the testing variable
to the model. Users specifiy the percent (not percentage points!) of improvement in R² (or the correlation) that a
squared term of a testing variable needs to add to be selected as a nonlinear auxiliary. For example, when this value is set
to 0.05 (5%), all squared terms of continuous predictors are considered that improve the R² by at least 5%.

{phang}{opt catlimit(#)} Users who need to check a lot of candidates as auxiliaries potentially might not have the time to declare
every single categorical one with the i.prefix. This option serves as a heuristic for the detection of categorical testing variables.
If the number of levels is equal to or below the set limit (the default value is 7), this testing variable will be declared with the i.prefix
and handled as a categorical one. Of course, this rule of thumb has severe limitations (e.g. the declaration of geographical
clusters, which can have hundreds of levels). Users should be aware of this and declare the relevant variables manually.

{phang}{opt saving(filename)} saves the final results as a .dta file.

{phang}{opt folds(#)} is the number of folds for the lasso cross validation procedure (CV). The higher the number of folds, the more stable
and precise the findings. While larger values are better, this increases computational time. For details refer to {helpb lasso}.

{phang}{opt seed(#)} sets the random-number seed. Because lasso relies on randomness, results may vary even with identical data
over multiple calls to {cmd:auxfinder}. Setting a seed ensures reproducibility.

{phang}{opt skiplasso} omits the computation of lasso models for final validation. This increases the speed of the computation manifold
as the lasso models are by far the longest running ones of the command. Use this option if the computational times are unacceptable.
Of course, then the final results are based solely on correlations and case numbers.

{phang}{opt elasticnet} computes elasticnet validation models instead of lasso ones. For details see {helpb elasticnet}.

{phang}{opt lassoopts(string)} is implemented so users can pass arbitrary options to customize the computation of
lasso / elasticnet models.

{phang}{opt noisily} shows more output, of the regression and lassomodels to give users additional insight of what is computed.

{title:Examples}

{dlgtab:Finding auxiliaries}

{p 8 12 2}. {stata "webuse nhanes2, clear"}{p_end}
{p 8 12 2}. {stata "replace hlthstat = . if hlthstat == 8"}{p_end}
{p 8 12 2}. {stata "auxfinder hlthstat albumin vitaminc copper hdresult, testvars(sampl - lead) seed(123)"}{p_end}


{dlgtab:Testing quadratic influences}

{p 8 12 2}. {stata "auxfinder hlthstat albumin vitaminc copper hdresult, testvars(sampl - lead) seed(123) nonlinear(0.05)"}{p_end}

{title:Returned results}

{pstd}Scalars:

{p2colset 5 27 20 2}{...}
{p2col : {cmd:r(n_strong_squared)}}number of strong squared auxils{p_end}
{p2col : {cmd:r(n_strong_cont}}number of strong continuous auxils{p_end}
{p2col : {cmd:r(n_strong_cat)}}number of strong categorical auxils{p_end}
{p2col : {cmd:r(n_weak_squared)}}number of weak squared auxils{p_end}
{p2col : {cmd:r(n_weak_cont)}}number of weak continuous auxils{p_end}
{p2col : {cmd:r(n_weak_cat)}}number of weak categorical auxils{p_end}
{p2col : {cmd:r(n_total_tested)}}number of variables passing the first check{p_end}


{pstd}Macros:

{p2col : {cmd:r(testvars)}}variables tested{p_end}
{p2col : {cmd:r(weak_squared)}}squared terms selected as weak auxils{p_end}
{p2col : {cmd:r(weak_cont)}}continuous variables selected as weak auxils{p_end}
{p2col : {cmd:r(weak_cat)}}categorical variables selected as weak auxils{p_end}
{p2col : {cmd:r(strong_squared)}}squared terms selected as strong auxils{p_end}
{p2col : {cmd:r(strong_cont)}}continuous variables selected as strong auxils{p_end}
{p2col : {cmd:r(strong_cat)}}categorical variables selected as strong auxils{p_end}

{title:References}

{phang}Bittmann, F. (2026). Auxfinder - Finding auxiliary variables for
imputation models {it:Zenodo.} https://doi.org/10.5281/zenodo.18536805

{title:Author}

{pstd}Felix Bittmann, Leibniz Institute for Educational Trajectories (LIfBi), felix.bittmann@lifbi.de

{pstd}Thanks for citing this software as follows:

{pmore}
Bittmann, F. (2026). {it:auxfinder: Finding auxiliary variables for imputation models.} https://doi.org/10.5281/zenodo.18536805

{title:Also see}

{helpb mi impute}, {helpb lasso}, {helpb elasticnet}
