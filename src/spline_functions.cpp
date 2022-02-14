#include <Rcpp.h>
using namespace Rcpp;

double cubic_product_int(double knot1, double knot2, double to, NumericVector coef1, NumericVector coef2) {
  if(knot1 > knot2) {
    NumericVector tmp_coef = clone(coef1);
    double tmp_knot = knot1;
    coef1 = clone(coef2);
    coef2 = clone(tmp_coef);
    knot1 = knot2;
    knot2 = tmp_knot;
  }
  double diff1 = to - knot1;
  double diff2 = to - knot2;

  double diff1_2 = pow(diff1, 2);
  double diff2_2 = pow(diff2, 2);
  double diff2_3 = pow(diff2, 3);
  double diff2_4 = pow(diff2, 4);
  double diff2_5 = pow(diff2, 5);
  double diff2_6 = pow(diff2, 6);

  double s = (coef1[0]+coef1[1]*diff1+coef1[2]*diff1_2+coef1[3]*pow(diff1, 3))*
    (coef2[0]*diff2+0.5*coef2[1]*diff2_2+1.0/3*coef2[2]*diff2_3+0.25*coef2[3]*diff2_4);
  s -= 0.5*(coef1[1]+2*coef1[2]*diff1+3*coef1[3]*diff1_2)*
    (coef2[0]*diff2_2+1.0/3*coef2[1]*diff2_3+1.0/6*coef2[2]*diff2_4+0.1*coef2[3]*diff2_5);
  s += 1.0/6*(2*coef1[2]+6*coef1[3]*diff1)*
    (coef2[0]*diff2_3+0.25*coef2[1]*diff2_4+0.1*coef2[2]*diff2_5+0.05*coef2[3]*diff2_6);
  s -= 1.0/24*6*coef1[3]*(coef2[0]*diff2_4+0.2*coef2[1]*diff2_5+1.0/15*coef2[2]*diff2_6+1.0/35*coef2[3]*pow(diff2, 7));
  return s;
}

IntegerVector find_interval(NumericVector v, NumericVector x) {
  Rcpp::IntegerVector res(v.length());
  for(int i=0; i < res.length(); ++i) {
    res[i] = std::distance(x.begin(), std::upper_bound(x.begin(), x.end(), v[i]));
  }
  return res;
}

double l2_inner_product(NumericVector knots_1, NumericMatrix coef_1, NumericVector knots_2, NumericMatrix coef_2, double from, double to) {
  NumericMatrix new_coef_1(coef_1.nrow() + 1, coef_1.ncol());
  NumericMatrix new_coef_2(coef_2.nrow() + 1, coef_2.ncol());

  new_coef_1(0, _) = NumericVector::create(coef_1(0,0)-knots_1(0)*coef_1(0, 1), coef_1(0, 1),0,0);
  new_coef_2(0, _) = NumericVector::create(coef_2(0,0)-knots_2(0)*coef_2(0, 1), coef_2(0, 1),0,0);
  for(int i=1; i<coef_1.nrow() + 1; i++) {
    new_coef_1(i, _) = coef_1(i-1, _);
  }
  for(int i=1; i<coef_2.nrow() + 1; i++) {
    new_coef_2(i, _) = coef_2(i-1, _);
  }

  coef_1 = new_coef_1;
  coef_2 = new_coef_2;

  NumericVector new_knots_1(knots_1.length()+2);
  NumericVector new_knots_2(knots_2.length()+2);

  new_knots_1[0] = from;
  new_knots_2[0] = from;
  new_knots_1[Range(1, knots_1.length())] = knots_1;
  new_knots_2[Range(1, knots_2.length())] = knots_2;
  new_knots_1[knots_1.length()+1] = to;
  new_knots_2[knots_2.length()+1] = to;

  knots_1 = new_knots_1;
  knots_2 = new_knots_2;

  NumericVector knots(knots_1.length()+knots_2.length());
  knots[Range(0, knots_1.length()-1)] = knots_1;
  knots[Range(knots_1.length(), knots_1.length()+knots_2.length()-1)] = knots_2;
  knots = sort_unique(knots);

  IntegerVector which_knot_1 = find_interval(knots, knots_1)-1;
  IntegerVector which_knot_2 = find_interval(knots, knots_2)-1;


  double inner_prod = 0;
  for(int i = 0; i < knots.length()-1; i++) {
    inner_prod += cubic_product_int(knots_1[which_knot_1[i]], knots_2[which_knot_2[i]],
                    knots[i+1], coef_1(which_knot_1[i], _), coef_2(which_knot_2[i], _));
  }

  return inner_prod;
}

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





