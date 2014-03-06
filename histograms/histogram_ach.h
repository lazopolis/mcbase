/** \file CHistogram.h
  *
  * ...
  */

#ifndef CHISTOGRAM_H
#define CHISTOGRAM_H
#include <iostream>
#include <sstream>
#include <fstream>
using namespace std;
#include "event.h"
#include "bin.h"
#include "user_interface.h"

class SimpleHistogram{
public:
    SimpleHistogram(unsigned, unsigned, const double&, const double &,
               const std::string &);
    ~SimpleHistogram(){};
    void bin_val(const double& x,const double &w);
    string give_name(){return _name;}
    void set_parameters(int NB,int LE, int HE);
    unsigned size() const { return _numbins; }
    void end(unsigned);
    void update();

    string info(){
        stringstream stream;stream<<give_name()
        <<" : "<<_numbins-1<<" bins\t["<<_lowend
        <<","<<_highend<<"]";return stream.str();
    }
    
    string plotinfo();
    void flush(){for (int i=0;i<_all_bins.size();i++) _all_bins[i].flush();}
    string xml();
    friend ostream& operator<<(ostream&, const SimpleHistogram&);
    friend string compare_histograms( const SimpleHistogram* H1,const SimpleHistogram* H2,const string& comp_type);
    friend string compare_histograms( const vector<SimpleHistogram*> &);
    friend string compare_histograms( const vector<SimpleHistogram*> &,bool color);

    double tot_running_f();
    unsigned    _numbins;
    unsigned    _firstbin;
    double      _lowend;
    double      _highend;
    std::string _name;
    double      _binsize;
    std::vector<Bin> _all_bins;

};



/** \brief Histogramming class
  *
  * ... */
class CHistogram: public SimpleHistogram
{
protected:
 

public:
  /** \brief Constructor 
   *
   * Initialize all values ......... */
    CHistogram(unsigned numbins_, unsigned firstbin_,
               const double& lowend_, const double& highend_,
               const std::string& name_, bool adapt_)
    : SimpleHistogram( numbins_, firstbin_,
                      lowend_,highend_,
                      name_){};
  /** \brief Empty destructor */
  ~CHistogram() {}

   /** \brief ... */
  void bin_event(const CombinedEvent&,const double &);
  /** \brief */
  virtual double determine_xval(const CombinedEvent&) = 0;
    
    friend string compare_histograms( const CHistogram* H1,const CHistogram* H2,const string& comp_type);
    friend string compare_histograms( const vector<CHistogram*> &);
    friend string compare_histograms( const vector<CHistogram*> &,bool color);
    
    bool the_event_is_in_ith_bin(int bin_id,const CombinedEvent&);
  void set_bin(int i,const double&res,const double& err,const double& prob){_all_bins[i].set(res,err,prob);}
};

#include "user_defined_histograms.h"




class HistogramBox
{
public://methods
    HistogramBox(const UserInterface & UI);
    void show_histogram_info_and_exit();
    void  book_histograms(const CombinedEvent&,const double & vegas_weight);
    void update_histograms_end_of_iteration(int NOP);
	void update_histograms_end_of_vegas_point();
    void print_histograms();
    string print_histograms_to_string();
    void write_to_histogram_file();
    CHistogram* ptr_to_histogram_with_id(unsigned m){return histogram_vector[m];}
    int size(){return histogram_vector.size();}
    
private://data
    vector<CHistogram*> histogram_vector;
    vector<CHistogram*> available_histograms;
    
    
private://methods
};




#endif
