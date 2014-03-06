/** \file CBin.cpp
  *
  * Implement CBin.h
  */

#include <math.h>
#include <iostream>
#include <iomanip>
#include <sstream>
using namespace std;

#include "bin.h"



void Bin::update(unsigned NOP)
{
  if (running_f!=0.0)
  {
    // Cumulated sums...
    // fabs takes care of the situation in which the integrand is constant and NOP*f2-f^2 is marginally negative due to round-off
    const double sigma=sqrt(fabs((NOP*running_f2-running_f*running_f)/(NOP-1.0)));
    // If first iteration
    if (avg_err==0.0)
    {
      // Values after first iteration
      avg_f           = running_f;
      avg_err         = sigma;
      running_chi_sq  = 0.0;
      iteration_number= 1;
    }
    else
    {
      // Save values
      const double old_err=avg_err;
      const double old_f=avg_f;
      // Update
      avg_err            = 1.0/sqrt(1.0/(avg_err*avg_err) + 1.0/(sigma*sigma));
      avg_f              = (avg_f/(old_err*old_err) + running_f/(sigma*sigma))*avg_err*avg_err;

     ++iteration_number;
     running_chi_sq    += pow(running_f - avg_f, 2.0)/avg_err;
     chi_sq             = running_chi_sq*iteration_number/(iteration_number-1.0);
    }
    // Reset running variables
    running_f=0.0;
    running_f2=0.0;
  }
  else
     {
       ++iteration_number;
       running_chi_sq    += pow( avg_f, 2.0)/avg_err;
       chi_sq             = running_chi_sq*iteration_number/(iteration_number-1.0);
     }

}
/*
double Bin::give_chi_sq() const
{   
    const double x=chi_sq;
    const int k=iteration_number;
    return(pow(x,k/2-1.0)*exp(-x/2.0)/pow(2.0,k/2.0)/my_factorial(k/2.0));
}
*/

void Bin::end()
{
  // Add values
  running_f     += intra_point_f;
  running_f2    += pow(intra_point_f, 2.0);
  // Reset values
  intra_point_f  = 0.0;
}

void Bin::add_single_point_package(const double &w)
{
    intra_point_f += w;
    point_counter++;
    end();
    points_since_last_update++;
}

void Bin::end_of_iteration_update()
{
    update(points_since_last_update);
    points_since_last_update=0.0;
}


string Bin::xml()
{
    stringstream stream(stringstream::out);
    stream<<"\n\t\t\t<s>"<< avg_f<<"</s>";
    stream<<"\n\t\t\t<e>"<<avg_err<<"</e>";
    stream<<"\n\t\t\t<epc>"<< fabs(avg_err/avg_f)*100.0
            <<"</epc>";
    stream<<"\n\t\t\t<chi>"<< chi_sq<<"</chi>";
    stream<<"\n\t\t\t<dof>"<< iteration_number<<"</dof>";
    stream<<"\n\t\t\t<nop>"<<point_counter <<"</nop>";
    return stream.str();
}

ostream& operator<<(ostream& stream, const Bin& bin)
{
  
  return stream <<setw(12)<< bin.avg_f
                << " | " <<setprecision(2)<<setw(8)<< bin.avg_err
                <<setprecision(4)
                << " | "<<setw(8)<<fabs(bin.avg_err/bin.avg_f)*100.0<<" %"
                << " | " <<setw(12)<< bin.chi_sq
                <<" / "<<setw(3)<<bin.iteration_number
                <<" | "<<setw(6)<<bin.point_counter<<" points";
                ;
}









