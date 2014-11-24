#ifdef _OPENMP
#include <omp.h>
#endif
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include "cpploss.h"
#include <Rcpp.h>
#include <math.h>
#include <Rmath.h>
#include <iostream>  

// [[Rcpp::export]]
SEXP  GCPM_cpploss(SEXP default_distr_a,SEXP link_function_a, SEXP S_a,SEXP Sigma_a, SEXP W_a, SEXP PD_a, SEXP PL_a, SEXP calc_rc_a, SEXP loss_thr_a,SEXP seed_a,SEXP max_entries_a){
  NumericMatrix S(S_a), W(W_a),Sigma(Sigma_a);
  NumericVector PD(PD_a), PL(PL_a),seed(seed_a),max_entries(max_entries_a),default_distr(default_distr_a),link_function(link_function_a),calc_rc(calc_rc_a),loss_thr(loss_thr_a);

  bool full=false;
  int N=S.nrow(), NS=S.ncol();
  int NCP=W.nrow(),defaults=0;
  int print_every=(int)floor(N/100),VaRszen=0,MaxVaRszen=std::min((int)(floor(max_entries(0)/NCP)),1000),inflate=1000,temp_int_1=0, temp_int_2=0;
  double condPD=0,temp=0,temp2=0;
  
  if(calc_rc(0)==1){
    temp_int_1=NCP;
    temp_int_2=MaxVaRszen;
  }
  else{
    temp_int_1=1;
    temp_int_2=1;
  }
  double *lossszenarios=new double[temp_int_2];
  double *buffer=new double[temp_int_2];
  double **buffer2=new double *[temp_int_1];
  double **CPsimlosses=new double *[temp_int_1];
  for(int i=0;i<temp_int_1;i++){
    CPsimlosses[i]=new double[temp_int_2];
    buffer2[i]=new double[temp_int_2];
  }
      
  double *CPlosses=new double[NCP];
  List ret;
  NumericVector simlosses(N),state(1);
  double *DD=new double[NCP];
  if(link_function(0)==2){
    for(int i=0;i<NCP;i++)
      DD[i]=R::qnorm(PD[i],0,1,1,0);
  }
  
  GetRNGstate();
  Progress p((int)floor(N/print_every), true);  
  if(link_function(0)==1){ //CRP
    for(int n=0;n<N;n++){
      if(n%print_every==0 && n>0){
        p.increment();
      }
      if(!Progress::check_abort()){
        for(int i=0;i<NCP;i++){
          condPD=W(i,0)*PD[i];
          for(int k=0;k<NS;k++)
            condPD+=W(i,k+1)*S(n,k)*PD[i];
          if(default_distr(i)==1){ //Bernoulli
            condPD=std::min(1.0,std::max(0.0000001,condPD));
            temp=unif_rand();
            if(temp<=condPD){
              CPlosses[i]=PL(i);
              simlosses(n)+=PL(i);
            }
            else
              CPlosses[i]=0;
          }
          else{ //Poisson
            condPD=std::max(0.0000001,condPD);
            
            defaults=(int)R::rpois(condPD);
            simlosses(n)+=PL(i)*defaults;
            CPlosses[i]=PL(i)*defaults;
          }
        }
        if(calc_rc(0)==1 && !full){
          if(simlosses(n)>=loss_thr(0)){
            VaRszen++;
            if(VaRszen>MaxVaRszen){
              if((MaxVaRszen+inflate)*NCP<=max_entries(0)){
                buffer=new double[MaxVaRszen];
                buffer2=new double *[NCP];
                for(int p=0;p<NCP;p++)
                  buffer2[p]=new double[MaxVaRszen];
                for(int s=0;s<MaxVaRszen;s++){
                  buffer[s]=lossszenarios[s];
                  for(int p=0;p<NCP;p++)
                    buffer2[p][s]=CPsimlosses[p][s];
                }
                MaxVaRszen+=inflate;
                lossszenarios=new double[MaxVaRszen];
                CPsimlosses=new double *[NCP];
                for(int i=0;i<NCP;i++)
                  CPsimlosses[i]=new double[MaxVaRszen];
                for(int s=0;s<(VaRszen-1);s++){
                  lossszenarios[s]=buffer[s];
                  for(int p=0;p<NCP;p++)
                    CPsimlosses[p][s]=buffer2[p][s];
                }
              }
              else{
                full=true;
              }
            }
            if(VaRszen<=MaxVaRszen){
              lossszenarios[VaRszen-1]=n;
              for(int i=0;i<NCP;i++)
                CPsimlosses[i][VaRszen-1]=CPlosses[i];
            }
            else{
              VaRszen--;
            }
          }
        }
      }
      else{
        simlosses(n)=-1;
      }
    }
  }
  else if(link_function(0)==2){ //CM
    for(int n=0;n<N;n++){
      if(n%print_every==0 && n>0){
        p.increment();
      }
      if(!Progress::check_abort()){
        for(int i=0;i<NCP;i++){
          temp=DD[i];
          for(int k=0;k<NS;k++){
            temp-=W(i,k)*S(n,k);
          }
          temp2=0;
          for(int k=0;k<NS;k++){
            for(int l=0;l<NS;l++)
              temp2+=W(i,k)*Sigma(k,l)*W(i,l);
          }
          condPD=R::pnorm(temp/sqrt(1-temp2),0,1,1,0);
          if(default_distr(i)==1){ //Bernoulli
            condPD=std::min(1.0,std::max(0.0000001,condPD));
            temp=unif_rand();
            if(temp<=condPD){
              CPlosses[i]=PL(i);
              simlosses(n)+=PL(i);
            }
            else
              CPlosses[i]=0;
          }
          else{ //Poisson
            condPD=std::max(0.0000001,condPD);
            defaults=(int)R::rpois(condPD);
            simlosses(n)+=PL(i)*defaults;
            CPlosses[i]=PL(i)*defaults; 
          }
        }
        if(calc_rc(0)==1  && !full){
          if(simlosses(n)>=loss_thr(0)){
            VaRszen++;
            if(VaRszen>MaxVaRszen){
              if((MaxVaRszen+inflate)*NCP<=max_entries(0)){
                buffer=new double[MaxVaRszen];
                buffer2=new double *[NCP];
                for(int p=0;p<NCP;p++)
                  buffer2[p]=new double[MaxVaRszen];
                for(int s=0;s<MaxVaRszen;s++){
                  buffer[s]=lossszenarios[s];
                  for(int p=0;p<NCP;p++)
                    buffer2[p][s]=CPsimlosses[p][s];
                }
                MaxVaRszen+=inflate;
                lossszenarios=new double[MaxVaRszen];
                CPsimlosses=new double *[NCP];
                for(int i=0;i<NCP;i++)
                  CPsimlosses[i]=new double[MaxVaRszen];
                for(int s=0;s<(VaRszen-1);s++){
                  lossszenarios[s]=buffer[s];
                  for(int p=0;p<NCP;p++)
                    CPsimlosses[p][s]=buffer2[p][s];
                }
              }
              else{
                full=true;
              }
            }
            if(VaRszen<=MaxVaRszen){
              lossszenarios[VaRszen-1]=n;
              for(int i=0;i<NCP;i++)
                CPsimlosses[i][VaRszen-1]=CPlosses[i];
            }
            else{
              VaRszen--;
            }
          }
        }
      }
      else{
        simlosses(n)=-1;
      }
    }
  }
  PutRNGstate();
  p.increment();

  NumericMatrix CPsimlossesfinal(NCP,VaRszen);
  NumericVector lossszenariosfinal(VaRszen);
  if(calc_rc(0)==1){
    if(full){
      for(int i=0;i<VaRszen;i++){
        for(int j=0;j<NCP;j++)
          CPsimlossesfinal(j,i)=-1;
        lossszenariosfinal(i)=-1;
      }
    }
    else{
      for(int i=0;i<VaRszen;i++){
        for(int j=0;j<NCP;j++)
          CPsimlossesfinal(j,i)=CPsimlosses[j][i];
        lossszenariosfinal(i)=lossszenarios[i]+1;
      } 
    }  
  }
  ret["simlosses"]=simlosses;  
  ret["CPsimlosses"]=CPsimlossesfinal;
  ret["lossszenarios"]=lossszenariosfinal;
  
  delete[] lossszenarios;
  delete[] buffer;
  delete[] CPlosses;
  for(int i=0;i<temp_int_1;i++){
    delete[] CPsimlosses[i];
    delete[] buffer2[i];
  }
  delete[] CPsimlosses;
  delete[] buffer2;
  return ret;
}

