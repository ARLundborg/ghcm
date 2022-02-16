#include <Rcpp.h>
using namespace Rcpp;
#include <algorithm>



//' Letting knot1 < knot2, this function computes \int_knot2^to p1(x-knot1) p2(x-knot2) dx
//' where p1 and p2 are third order polynomials with coefficients coef1 and coef2 respectively.
//' If knot2 > knot1 the labels are switched.
//'
//' @param knot1,knot2 shift values of the polynomials.
//' @param to upper limit of the integral.
//' @param coef1,coef2 coefficients of the third order polynomials.
//' @return the integral as described above.
double cubic_product_int(double knot1, double knot2, double to,
                             double coef1_0, double coef1_1, double coef1_2, double coef1_3,
                             double coef2_0, double coef2_1, double coef2_2, double coef2_3) {
  if(knot1 > knot2) {
    std::swap(knot1, knot2);
    std::swap(coef1_0, coef2_0);
    std::swap(coef1_1, coef2_1);
    std::swap(coef1_2, coef2_2);
    std::swap(coef1_3, coef2_3);
  }
  double diff1 = to - knot1;
  double diff2 = to - knot2;

  double diff1_2 = diff1*diff1;
  double diff1_3 = diff1_2*diff1;
  double diff2_2 = diff2*diff2;
  double diff2_3 = diff2_2*diff2;
  double diff2_4 = diff2_3*diff2;
  double diff2_5 = diff2_4*diff2;
  double diff2_6 = diff2_5*diff2;
  double diff2_7 = diff2_6*diff2;

  double one_third = 1.0/3;

  double s = (coef1_0+coef1_1*diff1+coef1_2*diff1_2+coef1_3*diff1_3)*
    (coef2_0*diff2+0.5*coef2_1*diff2_2+one_third*coef2_2*diff2_3+0.25*coef2_3*diff2_4);
  s -= 0.5*(coef1_1+2*coef1_2*diff1+3*coef1_3*diff1_2)*
    (coef2_0*diff2_2+one_third*coef2_1*diff2_3+0.5*one_third*coef2_2*diff2_4+0.1*coef2_3*diff2_5);
  s += 0.5*one_third*(2*coef1_2+6*coef1_3*diff1)*
    (coef2_0*diff2_3+0.25*coef2_1*diff2_4+0.1*coef2_2*diff2_5+0.05*coef2_3*diff2_6);
  s -= 0.25*coef1_3*(coef2_0*diff2_4+0.2*coef2_1*diff2_5+0.2*one_third*coef2_2*diff2_6+1.0/35*coef2_3*diff2_7);
  return s;
}

//' Computes the integral of the product of two natural cubic splines with knot
//'  sequences knots_1 and knots_2 and coefficient matrices coef_1 and coef_2
//'  from the splines::interpSpline function in R. The splines are assumed to live
//'  in a function space on the interval [from, to].
//'
//' @param knots_1,knots_2 knot sequences for the natural cubic splines.
//' @param coef_1,coef_2 coefficient matrices for the natural cubic splines.
//' @param from,to limits of integration.
//' @return the inner product of the splines.
double l2_inner_product(NumericVector knots_1, NumericMatrix coef_1, NumericVector knots_2, NumericMatrix coef_2, double from, double to) {

  double inner_prod = 0;
  int i=-1;
  int j=-1;
  double knot1 = from;
  double knot2 = from;
  double next_knot1 = knots_1[0];
  double next_knot2 = knots_2[0];
  double upper_lim = from;
  double coef1_0, coef1_1, coef1_2, coef1_3, coef2_0, coef2_1, coef2_2, coef2_3;
  while(upper_lim < to) {
    if(next_knot1 < next_knot2) {
      upper_lim = next_knot1;
    } else {
      upper_lim = next_knot2;
    }

    if (i > -1) {
      coef1_0 = coef_1(i, 0);
      coef1_1 = coef_1(i, 1);
      coef1_2 = coef_1(i, 2);
      coef1_3 = coef_1(i, 3);
    } else {
      coef1_0 = coef_1(0, 0) - coef_1(0, 1)*(knots_1[0] - from);
      coef1_1 = coef_1(0, 1);
      coef1_2 = 0;
      coef1_3 = 0;
    }

    if (j > -1) {
      coef2_0 = coef_2(j, 0);
      coef2_1 = coef_2(j, 1);
      coef2_2 = coef_2(j, 2);
      coef2_3 = coef_2(j, 3);
    } else {
      coef2_0 = coef_2(0, 0) - coef_2(0, 1)*(knots_2[0] - from);
      coef2_1 = coef_2(0, 1);
      coef2_2 = 0;
      coef2_3 = 0;
    }

    inner_prod += cubic_product_int(knot1, knot2, upper_lim,
                                    coef1_0, coef1_1, coef1_2, coef1_3,
                                    coef2_0, coef2_1, coef2_2, coef2_3);

    if(next_knot1 < next_knot2) {
      i++;
      knot1 = next_knot1;
      if (i < knots_1.length()-1) {
        next_knot1 = knots_1[i+1];
      } else {
        next_knot1 = to;
      }
    } else {
      j++;
      knot2 = next_knot2;
      if (j < knots_2.length()-1) {
        next_knot2 = knots_2[j+1];
      } else {
        next_knot2 = to;
      }
    }
  }

  return inner_prod;
}

//' Computes the matrix of L2 inner products of the splines given in list_of_splines
//' as produced by splines::interpSpline. The splines are assumed to be
//' functions on the interval [from, to].
//'
//' @param list_of_splines list of interpSpline objects.
//' @param from,to limits of integration.
//' @return matrix of inner products.
// [[Rcpp::export]]
NumericMatrix inner_product_matrix_splines(List list_of_splines, double from, double to) {
  int n = list_of_splines.size();
  NumericMatrix mat(n, n);
  for(int i=0; i<n; i++) {
    for(int j=0; j<n; j++) {
      if(i > j) {
        mat(i, j) = mat(j, i);
      } else {
        List spline_i = list_of_splines[i];
        List spline_j = list_of_splines[j];
        NumericVector knots_i = spline_i["knots"];
        NumericVector knots_j = spline_j["knots"];
        NumericMatrix coef_i = spline_i["coefficients"];
        NumericMatrix coef_j = spline_j["coefficients"];
        mat(i, j) = l2_inner_product(knots_i, coef_i, knots_j, coef_j, from, to);
      }
    }
  }
  return mat;
}
