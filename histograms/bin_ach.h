/** \file CBin.h
  *
  * Declares Bin object used in the CHistogram class.
  */

#ifndef BIN_H
#define BIN_H



/** \brief Bin of histogram
  *
  * Represent a bin of a monte-carlo histogram. The bins are in fact integrating themselves since the VEGAS adaptation is in general only done for the total cross section. The weight thus have to passed from the VEGAS generator.
  * 
  * A single VEGAS point must usually be split in different pieces (usually corresponding to different subtraction terms), and thus the computation of errors, which involves the square of the value at each point, has to performed after all pieces of that point have been added with the add() function. When all the pieces of a point have been added the function end() should be called to update the values and its squared. */
class Bin {
    /** \brief Accumulated value */
    double running_f;
    /** \brief Accumulated square */
    double running_f2;
    /** \brief Accumulated chi-squared */
    double running_chi_sq;
    /** \brief Vegas average */
    double avg_f;
    /** \brief Average error (stand. dev.) */
    double avg_err;
    /** \brief Chi-squared*/
    double chi_sq;
    /** \brief */
    double intra_point_f;
    /** \brief Counter */
    unsigned iteration_number;
    //* \brief number of points in this bin */
    unsigned point_counter;
  public:
    void flush(){running_f=0.0;running_f2=0.0;running_chi_sq=0.0;avg_f=0.0;avg_err=0.0;chi_sq=0.0;intra_point_f=0.0;iteration_number=0;point_counter=0;}
    string xml();
    void set(const double& res,const double& err,const double& prob){avg_f = res;avg_err = err;chi_sq = prob;}
    double give_chi_sq() const;
    /** \brief Constructor
    *
    * Puts all internal variables to 0. */
  Bin()
    : running_f(0.0), running_f2(0.0), avg_f(0.0), avg_err(0.0), intra_point_f(0.0), running_chi_sq(0.0), iteration_number(0),point_counter(0),
            points_since_last_update(0)
  {}

  /** \brief Accumulate point */
  void add(const double& w)
    { intra_point_f += w;point_counter++; }
  /** \brief Update all internals
    *
    * Computes the new value error and chi-squared of the integral. 
    * Parameter is the number of points that have been thrown by VEGAS 
    * since tha last call to update(). */
  void update(unsigned);
  /** \brief End of point accumulation
    *
    * This functions ends and restarts the accumulation of points
    * belonging to the same physical point. */
  void end();
    
    
    void add_single_point_package(const double &w);
    void end_of_iteration_update();
    
  /** \brief Return value */
  operator double() const
  { return avg_f; }
  /** \brief Return value */
  double value() const
  { return avg_f; }
  /** \brief Return error */
  double error() const
  { return avg_err; }

  friend ostream& operator<<(ostream&, const Bin&);
    
  double give_running_f(){return running_f;}  
private://data
    unsigned points_since_last_update;
};

#endif
