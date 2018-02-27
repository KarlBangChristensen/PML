/*************************************************************************
macro to compute Martin-Löf test statistic for test of unidimensionality

	%pml(data, scale1, scale2, max, out=PML, bootstrap=NO, nsimu=1000);

the data set 'data' contains the items in 'scale1' and the items in 
'scale2' all items are scored 0-'max'. The macro can compute p-value 
suing parametric bootstrap
*************************************************************************/

%macro pml(data, scale1, scale2, max, out=PML, bootstrap=NO, nsimu=1000);

options nonotes nostimer; 
ods listing close;

*options mprint;
/* help macro to fit Rasch model using CML */
%macro r(data, items, max, parfile=_pf, gammafile=_gf, loglfile=_lf);
 options nonotes nostimer;
 ods listing close; 
/* save number of respondents and item names as macro variables */
 data _null_; 
  set &data end=final; 
  if final then call symput('N',trim(left(_N_))); 
 run;

 data _null_; 
  set &data;
  array _y (*) &items; 
  length name $8;
  if _n_=1 then do;
   do _i=1 to dim(_y);
    call vname(_y{_i},name); 
    call symput('_item'||trim(left(put(_i,4.))),trim(left(name)));
   end;
   _p=dim(_y); 
   call symput('_nitems',trim(left(put(_p,4.))));
  end;
 run;

/* constants */
 %let maxplus1=%eval(&max+1);
 %let maxscore=%eval(&max*&_nitems);



/* check score distribution */

%do _sc=0 %to &maxscore;
 %let _n&_sc=0;
%end;

data _null_; 
 set &data;
 t=0 %do _i=1 %to &_nitems; +&&_item&_i %end;;
 %do _sc=0 %to &maxscore;
  if t=&_sc then call symput("_n&_sc",1);
 %end;
run;

 

/* make table - code dummy variables */
proc freq data=&data; 
 table %do _i=1 %to &_nitems-1; &&_item&_i* %end; &&_item&_nitems
  /sparse out=_table noprint; 
run;
data _table; 
 set _table; 
 t=0 %do _i=1 %to &_nitems; +&&_item&_i %end;;
 %do _it=1 %to &_nitems; 
  %do _cat=1 %to &max; 
   v&_it&_cat=(&&_item&_it=&_cat); 
  %end; 
 %end;
 %do _sc=0 %to &maxscore;
  n&_sc=(t=&_sc);
 %end;
run;

/* poisson model */

ods output genmod.parameterestimates=&parfile;
ods output Genmod.ConvergenceStatus=_conv;

proc genmod data=_table; 
  model count=%do _it=1 %to &_nitems; 
   %do _cat=1 %to &max; 
    v&_it&_cat 
   %end; 
  %end;
  %do _sc=1 %to %eval(&maxscore-1);
   %if &&_n&_sc=1 %then %do; n&_sc %end;
  %end;
  /d=p noint converge=0.000000000001 convh=1 link=log maxiter=20000 itprint; 
 run;



/* standardize item parameter estimates, create data set, save as macro variables */
 data &parfile; 
  set &parfile; 
  keep parameter estimate; 
  if (parameter notin ('Intercept','t','Scale')); 
 run;
 
 proc transpose data=&parfile out=&parfile; 
  id parameter; 
 run;
 
 data &parfile; 
  set &parfile; 
  s=(%do _it=1 %to &_nitems-1; v&_it&max+ %end; v&_nitems&max)/(&_nitems*&max);
  %do _it=1 %to &_nitems; 
   %do _cat=1 %to &max; 
    eta&_it&_cat=v&_it&_cat-&_cat*s; 
   %end; 
  %end;
  keep %do _it=1 %to &_nitems; %do _cat=1 %to &max; eta&_it&_cat %end; %end;;
  %do _it=1 %to &_nitems; 
   call symput("eta&_it.0",trim(left(eta&_it|0))); 
  %end; 
  %do _it=1 %to &_nitems; 
   %do _cat=1 %to &max; 
    call symput("eta&_it.&_cat",trim(left(eta&_it&_cat)));
   %end; 
  %end;
 run;

/* compute gammas values by defining temporary gammas starting with item1 and
compute recursively for (item1,item2), (item1,item2,item3), ... save as macro 
variables and in file */
 data _gamma; 
  set &parfile;
  %do _it=1 %to &_nitems; 
   array g&_it._[0:%eval(&_it*&max)]; 
  %end;
  g1_[0]=1; 
  %do _c=1 %to &max; 
   g1_[&_c]=exp(&&eta1&_c); 
  %end; 
  %do _m=2 %to &_nitems; 
   do c&_m=0 to &_m*&max; 
    g&_m._[c&_m]=0; 
    %do _cat=0 %to &max;
     if &_cat<=c&_m<=%eval(&_m-1)*&max+&_cat then 
     g&_m._[c&_m]=g&_m._[c&_m]+exp(&&eta&_m&_cat)*g%eval(&_m-1)_[c&_m-&_cat];
    %end;
   end;
  %end;
  %do _sc=0 %to (&_nitems*&max); 
   call symput("gam&_sc",trim(left(g&_nitems._[&_sc]))); 
  %end;
 run;

 data &gammafile; 
  %do _sc=0 %to (&_nitems*&max); 
   gamma&_sc=&&gam&_sc; 
   lgamma&_sc=log(&&gam&_sc); 
  %end; 
 run;

/* compute maximum of log likelihood - save value as macro variable */
 data &loglfile;
  set &data;
  %do _it=1 %to &_nitems; 
   %do _cat=1 %to &max; 
    v&_it&_cat=(&&_item&_it=&_cat); 
   %end; 
  %end;
  sc=%do _i=1 %to &_nitems-1; &&_item&_i+ %end; &&_item&_nitems;
  pr=exp(0 %do _it=1 %to &_nitems; 
   %do _cat=1 %to &max; 
    +v&_it&_cat*&&eta&_it.&_cat 
   %end; 
  %end;)/(0 %do _sc=0 %to (&_nitems*&max); +&&gam&_sc*(sc=&_sc) %end;);
  lp=log(pr); 
  keep lp; 
 run;

 proc transpose data=&loglfile out=&loglfile; 
 run;

 data &loglfile; 
  set &loglfile;
  array l[&N] col1-col&N; 
  logl=0; 
  do i=1 to &N; 
   logl=logl+l[i]; 
  end; 
  keep logl; 
  call symput('ll',trim(left(logl)));
 run;

 /* end of help macro */
%mend r;





/* save number of respondents as macro variable */
data _null_; 
 set &data end=final; 
 if final then call symput('_N',trim(left(_N_))); 
run;

/* save the item names for first dimension as macro variables */
data _null_; 
 set &data; 
 array _y (*) &scale1; 
 length name $8;
 if _n_=1 then do;
  do _i=1 to dim(_y);
   call vname(_y{_i},name); 
   call symput('_scale1'||trim(left(put(_i,4.))),trim(left(name)));
  end;
  _p=dim(_y); 
  call symput('_nitems1',trim(left(put(_p,4.))));
 end; 
run;

/* save the item names for second dimension as macro variables */
data _null_; 
 set &data; 
 array _y (*) &scale2; 
 length name $8;
 if _n_=1 then do;
  do _i=1 to dim(_y);
   call vname(_y{_i},name); 
   call symput('_scale2'||trim(left(put(_i,4.))),trim(left(name)));
  end;
  _p=dim(_y); 
  call symput('_nitems2',trim(left(put(_p,4.))));
 end; 
run;

/* constants */
%let _maxsc1=%eval(&_nitems1*&max); 
%let _maxsc2=%eval(&_nitems2*&max); 
%let _ncomb=%eval(%eval(&_maxsc1+1)*%eval(&_maxsc2+1));
%let _nscores=%eval(&_maxsc1+&_maxsc2+1);
%let _maxsc=%eval(&_maxsc1+&_maxsc2);

/* new data with sub scores and total scores */
data &data; 
 set &data; 
 _subsc1=0 %do _i=1 %to &_nitems1; +&&_scale1&_i %end;;
 _subsc2=0 %do _i=1 %to &_nitems2; +&&_scale2&_i %end;;
 _score=_subsc1+_subsc2; 
run;

/* make tables */
ods output crosstabfreqs=_table onewayfreqs=_dist;
proc freq data=&data; 
 tables _subsc1*_subsc2 _score/nocol norow nopercent nocum; 
run;

/* compute sums */
data _table; 
 set _table (where=(_type_='11')); 
 keep frequency; 
run;

data _dist; 
 set _dist; 
 keep frequency; 
run;

proc transpose data=_table out=_table; run;

data _table; 
 set _table; 
 array col[&_ncomb] col1-col&_ncomb; 
 _sum=0; 
 do i=1 to &_ncomb; 
  if col[i]>0 then _sum=_sum+col[i]*log(col[i]/&_N); 
 end; 
run;

proc transpose data=_dist out=_dist; run;

/* save marginal score distribution as macro variable */
data _null_; 
 set _dist;
 %do _sc=1 %to &_maxsc+1;
  call symput('_number'||trim(left(put(%eval(&_sc-1),4.))),trim(left(col&_sc)));
 %end;
run;

data _dist; 
 set _dist; 
 array col[&_nscores] col1-col&_nscores; 
 _sum=0; 
 do i=1 to &_nscores; 
  if col[i]>0 then _sum=_sum+col[i]*log(col[i]/&_N); 
 end; 
run;

data _new; 
 set _table _dist; 
 keep _sum; 
run; 

proc transpose data=_new out=_new; run;

data _new1; 
 set _new; 
 _dif=col2-col1; 
 keep _dif; 
run; 

/* fit rasch models using help macro */
%r(data=&data, items=&scale1, max=&max, parfile=_pf1, gammafile=_gf1, loglfile=_lf1);
%r(data=&data, items=&scale2, max=&max, parfile=_pf2, gammafile=_gf2, loglfile=_lf2);
%r(data=&data, items=&scale1 &scale2, max=&max, parfile=_pf, gammafile=_gf, loglfile=_lf);

data _new2; 
 set _lf _lf1 _lf2; 
run;

proc transpose data=_new2 out=_new2; run;

data _new2; 
 set _new2; 
 _logl=col2+col3-col1; 
 keep _logl; 
run;

data &out; 
 merge _new1 _new2; 
 pml=2*(_logl-_dif); 
 df=&_ncomb-&_nscores-1; 
 p=1-probchi(pml,df); 
 call symput('_pml',trim(left(pml))); 
 call symput('_df',trim(left(df))); 
 call symput('_p',trim(left(p))); 
run;

/* write to log window */
%put -------------------------------------------------------------------;
%put Testing unidimensionality reading data: &data;
%put subscales &scale1 (scores 0-&_maxsc1) and &scale2 (scores 0-&_maxsc2);
%put Martin-Löf test statistic is &_pml (df=&_df, p=&_p);
%put computing table of observed and expected counts ;
%put -------------------------------------------------------------------;


/* observed counts */
proc sort data=&data; 
 by _subsc2; 
run;

%do _sc2=0 %to &_maxsc2; 
 ods output freq.bygroup%eval(&_sc2+1)._subsc1.onewayfreqs=_obs&_sc2; 
%end;

proc freq data=&data; 
 table _subsc1/nocum nopercent; 
 by _subsc2; 
run;

data _zeros; 
 %do _sc=0 %to &_maxsc1; 
  _subsc1=&_sc; 
  output; 
 %end; 
 keep _subsc1; 
run; 

%do _sc2=0 %to &_maxsc2; 
 data _obs&_sc2; 
  merge _obs&_sc2 _zeros; 
  by _subsc1; 
  score1=_subsc1; 
  sc2_&_sc2=frequency; 
  if sc2_&_sc2=. then sc2_&_sc2=0; 
  keep score1 sc2_&_sc2; 
 run;
%end;

data _obs; 
 merge %do _sc2=0 %to &_maxsc2; _obs&_sc2 %end;; 
run;

 /* expected counts - gammas for total scale and subscale gammas */
data _g; 
 set _gf(drop=lgamma0-lgamma&_maxsc); 
 %do _sc=0 %to &_maxsc; 
  rename gamma&_sc=_g_&_sc; 
 %end; 
run;
 
data _gg; 
 set _pf;
 %do _it=1 %to %eval(&_nitems1+&_nitems2); 
  eta&_it.0=0; 
 %end; 
 %do _it=1 %to &_nitems1; 
  array g&_it._[0:%eval(&_it*&max)]; 
 %end;
 g1_[0]=1; 
 %do _c=1 %to &max; 
  g1_[&_c]=exp(eta1&_c); 
 %end; 
 %do _m=2 %to &_nitems1; 
  do c&_m=0 to &_m*&max; 
   g&_m._[c&_m]=0; 
    %do _cat=0 %to &max;
     if &_cat<=c&_m<=%eval(&_m-1)*&max+&_cat then 
      g&_m._[c&_m]=g&_m._[c&_m]+exp(eta&_m&_cat)*g%eval(&_m-1)_[c&_m-&_cat];
     %end; 
    end; 
   %end;
   %do _sc=0 %to &_maxsc1; 
    _g1_&_sc=g&_nitems1._[&_sc]; 
   %end; 
  %do _it=1 %to &_nitems2; 
   array gg&_it._[0:%eval(&_it*&max)]; 
  %end;
  gg1_[0]=1; 
  %do _c=1 %to &max; 
   gg1_[&_c]=exp(eta%eval(&_nitems1+1)&_c); 
  %end; 
  %do _m=2 %to %eval(&_nitems2); 
   do c&_m=0 to &_m*&max; 
    gg&_m._[c&_m]=0; 
   %do _cat=0 %to &max;
    if &_cat<=c&_m<=%eval(&_m-1)*&max+&_cat then 
    gg&_m._[c&_m]=gg&_m._[c&_m]
     +exp(eta%eval(&_nitems1+&_m)&_cat)*gg%eval(&_m-1)_[c&_m-&_cat];
   %end;
  end;
 %end;
 %do _sc=0 %to &_maxsc2; 
  _g2_&_sc=gg&_nitems2._[&_sc]; 
 %end; 
 keep %do _sc=0 %to &_maxsc1; _g1_&_sc %end; %do _sc=0 %to &_maxsc2; _g2_&_sc %end;; 
run;

data _exp; merge _g _gg _dist; 
 %do _sc1=0 %to &_maxsc1; 
  score1=&_sc1; 
  %do _sc2=0 %to &_maxsc2; 
   sc2_&_sc2=col%eval(&_sc1+&_sc2+1)*(_g1_&_sc1*_g2_&_sc2)/_g_%eval(&_sc1+&_sc2);
   prob_&_sc2=(_g1_&_sc1*_g2_&_sc2)/_g_%eval(&_sc1+&_sc2);
  %end; 
  output; 
 %end; 
 keep score1 sc2_0-sc2_&_maxsc2 prob_0-prob_&_maxsc2; 
run;


/* compute residuals */
data _o; 
 set _obs;
 %do _sc2=0 %to &_maxsc2;
  rename sc2_&_sc2=obs&_sc2;
 %end;
run;

data _e; 
 set _exp;
 %do _sc2=0 %to &_maxsc2;
  rename sc2_&_sc2=exp&_sc2;
 %end;
run;

data _res;  
 merge _o _e; 
 by score1;
 %do _sc=0 %to &_maxsc; 
  N&_sc=&&_number&_sc; 
 %end;
 array _N[&_nscores] N0-n&_maxsc;
 %do _sc2=0 %to &_maxsc2;
  totsc=score1+&_sc2;
  if score1+&_sc2 in (0,&_maxsc) then res&_sc2=0; else
   res&_sc2=(obs&_sc2-exp&_sc2)*(_N[totsc+1]*prob_&_sc2*(1-prob_&_sc2))**-.5;
 %end;
 keep score1 res0-res&_maxsc2;
run;


/* compute bootstrap p-values if required */
%if %upcase(%left(%trim(&bootstrap)))=YES %then %do;

%put using parametric bootstrap to compute p-value;
%put iteration history written to file &out._history;
%put -------------------------------------------------------------------;

%let seed=0;
 
/* score probabilities */
%do _sc=0 %to &_maxsc; 
 %if (&&_number&_sc=.) %then %let _number&_sc=0;
%end;
/* initiate sim results file */
data _simval;
 length simu simpml 8;
run; 

data _realgamma;
 set _gamma;
run;


/* start loop */
%let over=0;
%do _simu=1 %to &nsimu;

/* simulate data set under the null hypothesis of unidimensionality */
data _simu;
 set _realgamma;
 %do _it=1 %to %eval(&_nitems1+&_nitems2);
  eta&_it.0=0;
 %end;
 %do _h=0 %to &max;
  g1_%eval(&_h+1)=exp(eta1&_h);
 %end;
 do _j=1 to &_N;
  /* simulate score values */
  score=rantbl((&seed*&_simu) %do _sc=0 %to &_maxsc; ,&&_number&_sc/&_N %end;)-1;     
  /* simulate item values from conditional distribution */
  sc=score;
  %do _it=%eval(&_nitems1+&_nitems2) %to 2 %by -1; 
   array pgam&_it (*) %do _nr=1 %to %eval(&_it*&max+1); g&_it._&_nr %end;;
   array prgam&_it (*) %do _nr=1 %to %eval((&_it-1)*&max+1); g%eval(&_it-1)_&_nr %end;;
   /* probabilites */
   array _pprob&_it[%eval(&max+1)]; 
   %do _h=0 %to &max; 
    _pprob&_it[&_h+1]=0;
   %end;
   /* possible range */
   %do _h=0 %to &max;
    if ((sc-%eval(&_it-1)*&max le &_h) and (&_h le sc)) then 
    _pprob&_it[%eval(&_h+1)]=min(1,exp(eta&_it.&_h)*prgam&_it[sc-&_h+1]/pgam&_it[sc+1]);
   %end;
   _item&_it=rantbl((&seed*&_simu) %do _h=0 %to &max; ,_pprob&_it[&_h+1] %end;)-1;
   sc=sc-_item&_it;
  %end;
  _item1=sc;
  person=_j;
  output;  
 end;
 keep _item1-_item%eval(&_nitems1+&_nitems2);
run;

/* check data set */
%let _ok=1;
%let _error=' ';
%do _it=1 %to %eval(&_nitems1+&_nitems2);
 %do _c=0 %to &max;
  %let I_&_it._&_c=0;
 %end;
%end;

data _null_;
 set _simu;
 %do _it=1 %to %eval(&_nitems1+&_nitems2);
  %do _c=0 %to &max;
   if _item&_it=&_c then call symput("I_&_it._&_c",1);
  %end;
 %end;
run;

%do _it=1 %to %eval(&_nitems1+&_nitems2);
 %do _c=0 %to &max;
  %if &&I_&_it._&_c=0 %then %do;
   %let _error="_item &_it response &_c not used";
   %let _ok=0;
  %end;
 %end;
%end;

/* compute pml statistic for the simulated data set if simulated data is ok */
%if &_ok=1 %then %do;

/* sub scores and total scores */
data _simu; 
 set _simu; 
 _subsc1=0 %do _i=1 %to &_nitems1; +_item&_i %end;;
 _subsc2=0 %do _i=%eval(1+&_nitems1) %to %eval(&_nitems1+&_nitems2); +_item&_i %end;;
 _score=_subsc1+_subsc2; 
run;

/* tables */
ods output crosstabfreqs=_simtable onewayfreqs=_simdist;
proc freq data=_simu; 
 tables _subsc1*_subsc2 _score/nocol norow nopercent nocum; 
run;

/* sums */
data _simtable; 
 set _simtable (where=(_type_='11')); 
 keep frequency; 
run;

data _simdist; 
 set _simdist; 
 keep frequency; 
run;

proc transpose data=_simtable out=_simtable; run;

data _simtable; 
 set _simtable; 
 array col[&_ncomb] col1-col&_ncomb; 
 _sum=0; 
 do i=1 to &_ncomb; 
  if col[i]>0 then _sum=_sum+col[i]*log(col[i]/&_N); 
 end; 
run;

proc transpose data=_simdist out=_simdist; run;

data _simdist; 
 set _simdist; 
 array col[&_nscores] col1-col&_nscores; 
 _sum=0; 
 do i=1 to &_nscores; 
  if col[i]>0 then _sum=_sum+col[i]*log(col[i]/&_N); 
 end; 
run;

data _simnew; 
 set _simtable _simdist; 
 keep _sum; 
run; 

proc transpose data=_simnew out=_simnew; run;

data _simnew1; 
 set _simnew; 
 _dif=col2-col1; 
 keep _dif; 
run; 




/* fit rasch models using help macro */
%r(data=_simu, items=_item1-_item&_nitems1, max=&max, parfile=_pf1, gammafile=_gf1, loglfile=_simlf1);
data _null_;
 set _conv;
 if status>0 then do;
  call symput('_ok',0);
  call symput('_error','model for sub scale 1 did not converge');
 end;
run;
%r(data=_simu, items=_item%eval(1+&_nitems1)-_item%eval(&_nitems1+&_nitems2), max=&max, parfile=_pf2, gammafile=_gf2, loglfile=_simlf2);
data _null_;
 set _conv;
 if status>0 then do;
  call symput('_ok',0);
  call symput('_error','model for sub scale 2 did not converge');
 end;
run;
%r(data=_simu, items=_item1-_item%eval(&_nitems1+&_nitems2), max=&max, parfile=_pf, gammafile=_gf, loglfile=_simlf);
data _null_;
 set _conv;
 if status>0 then do;
  call symput('_ok',0);
  call symput('_error','model for total scale did not converge');
 end;
run;
data _simnew2; 
 set _simlf _simlf1 _simlf2; 
run;

%if &_ok=1 %then %do;
 proc transpose data=_simnew2 out=_simnew2; run;
 data _simnew2; 
  set _simnew2; 
  _logl=col2+col3-col1; 
  keep _logl; 
 run;
 data _simvalue; 
  merge _simnew1 _simnew2; 
  simpml=2*(_logl-_dif); 
  simu=&_simu;
  call symput('_simpml',left(simpml));
 run;
 proc append base=_simval new=_simvalue(keep=simu simpml);
 run;
 %put simulation number &_simu : simulated value &_simpml;
%end;

%end;

%if &_ok=0 %then %do; 
 %put simulation number &_simu : &_error;
 data _simvalue; 
  simpml=.; 
  simu=&_simu;
 run;
 proc append base=_simval new=_simvalue(keep=simu simpml);
 run;
%end;


%end;

/* create output files */
data _null_;
 set _simval end=last;
 retain over(0);
 over=over+(simpml>&_pml);
 number=_N_-1;
 if last then do;
  p_mc=over/number;
  call symput('_p_mc',left(p_mc));
 end;
run;

data &out; 
 set &out;
 p_boot=&_p_mc;
 nsimu=&nsimu;
 drop _dif _logl;
run;

data &out._history; 
 set _simval;
 if simu ne .;
run;

/* end of bootstrap part */
proc datasets nolist;
 delete _realgamma _simdist _simlf _simlf1 _simlf2 _simnew _simnew1 _simnew2 
_simtable _simu _simval _simvalue;
run;
quit;


%end;


/* end of macro */
ods listing;
title 'observed counts'; 
proc print data=_obs noobs; run;
title 'expected counts under the hypothesis of unidimensionality';
proc print data=_exp(drop=prob_0-prob_&_maxsc2) noobs; run;
title 'residuals';
proc print data=_res; run;

proc datasets nolist;
 delete _conv _dist _e _exp _g _gamma _gf _gf1 _gf2 _gg _lf _lf1 _lf2 _new _new1 _new2
_o _obs _obs0-_obs&_maxsc2 _pf _pf1 _pf2 _res _table _zeros;
run; quit;

options notes stimer; 
title ' ';
%mend pml;
