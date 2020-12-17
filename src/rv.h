#pragma once

#include <Rcpp.h>

class RV {
public:
  virtual Rcpp::S4 to_S4() const = 0;
};

class BetaRV : public RV {
private:
  double shape1_;    /** positive shape parameter */
  double shape2_;    /** positive shape parameter */
public:
  /**
   * Getter for the 1st shape parameter
   */
  double get_shape1() const { return shape1_; }

  /**
   * Getter for the 2nd shape parameter
   */
  double get_shape2() const { return shape2_; }

  /**
   * Setter for the 1st shape parameter
   * @param value new positive shape parameter value
   */
  void set_shape1(double value) {
    if (value <= 0) Rcpp::stop("BetaRV: invalid shape1 parameter");
    shape1_ = value;
  }

  /**
   * Setter for the 2nd shape parameter
   * @param value new positive shape parameter value
   */
  void set_shape2(double value) {
    if (value <= 0) Rcpp::stop("BetaRV: invalid shape2 parameter");
    shape2_ = value;
  }

  /**
   * Cast to an R S4 object
   * @return a BetaRV S4 object
   */
  Rcpp::S4 to_S4() const {
    Rcpp::S4 out("BetaRV");
    out.slot("shape1") = shape1_;
    out.slot("shape2") = shape2_;
    return out;
  }

  /**
   * Constructor from a BetaRV S4 object
   * @param obj a BetaRV S4 object
   */
  BetaRV(const Rcpp::S4 &obj) {
    if (!obj.is("BetaRV")) Rcpp::stop("BetaRV: invalid BetaRV S4 object");
    shape1_ = obj.slot("shape1");
    shape2_ = obj.slot("shape2");
  }

  /**
   * Constructor
   * @param alpha positive alpha shape parameter
   * @param beta positive beta shape parameter
   * @return a BetaRV instance
   */
  BetaRV(double shape1, double shape2) 
    : shape1_(shape1), shape2_(shape2) {}
};

class GammaRV : public RV {
private:
  double shape_;   /** positive shape parameter */
  double rate_;    /** positive rate parameter */
public:
  /**
   * Getter for shape parameter
   */
  double get_shape() const { return shape_; }

  /**
   * Getter for rate parameter
   */
  double get_rate() const { return rate_; }

  /**
   * Setter for shape parameter
   * @param value new positive shape parameter value
   */
  void set_shape(double value) {
    if (value <= 0) Rcpp::stop("GammaRV: invalid shape parameter");
    shape_ = value;
  }

  /**
   * Setter for rate parameter
   * @param value new positive rate parameter value
   */
  void set_rate(double value) {
    if (value <= 0) Rcpp::stop("GammaRV: invalid rate parameter");
    rate_ = value;
  }

  /**
   * Cast to an R S4 object
   * @return a GammaRV S4 object
   */
  Rcpp::S4 to_S4() const {
    Rcpp::S4 out("GammaRV");
    out.slot("shape") = shape_;
    out.slot("rate") = rate_;
    return out;
  }

  /**
   * Constructor from a GammaRV S4 object
   * @param obj a GammaRV S4 object
   */
  GammaRV(const Rcpp::S4 &obj) {
    if (!obj.is("GammaRV")) Rcpp::stop("GammaRV: invalid GammaRV S4 object");
    shape_ = obj.slot("shape");
    rate_ = obj.slot("rate");
  }

  /**
   * Constructor
   * @param shape positive shape parameter
   * @param rate positive rate parameter
   * @return a GammaRV instance
   */
  GammaRV(double shape, double rate) 
    : shape_(shape), rate_(rate) {}
};

class ShiftedNegBinomRV : public RV {
private:
  int size_;     /** positive size parameter */
  double prob_;  /** prob parameter */
public:
  /**
   * Getter for size parameter
   */
  double get_size() const { return size_; }

  /**
   * Getter for prob parameter
   */
  double get_prob() const { return prob_; }

  /**
   * Setter for size parameter
   * @param value new positive size parameter value
   */
  void set_size(int value) {
    if (value <= 0) Rcpp::stop("ShiftedNegBinomRV: invalid size parameter");
    size_ = value;
  }

  /**
   * Setter for prob parameter
   * @param value new prob parameter value
   */
  void set_prob(double value) {
    if (value <= 0 || value > 1) Rcpp::stop("ShiftedNegBinomRV: invalid prob parameter");
    prob_ = value;
  }
  
  /**
   * Cast to an R S4 object
   * @return a BetaRV S4 object
   */
  Rcpp::S4 to_S4() const {
    Rcpp::S4 out("BetaRV");
    out.slot("size") = size_;
    out.slot("prob") = prob_;
    return out;
  }

  /**
   * Constructor from a ShiftedNegBinomRV S4 object
   * @param obj a ShiftedNegBinomRV S4 object
   */
  ShiftedNegBinomRV(const Rcpp::S4 &obj) {
    if (!obj.is("ShiftedNegBinomRV")) Rcpp::stop("ShiftedNegBinomRV: invalid ShiftedNegBinomRV S4 object");
    size_ = obj.slot("size");
    prob_ = obj.slot("prob");
  }
  
  /**
   * Constructor
   * @param size positive shape parameter
   * @param rate positive rate parameter
   * @return a ShiftedNegBinomRV instance
   */
  ShiftedNegBinomRV(int size, double prob) 
    : size_(size), prob_(prob) {}
};