//
// Created by Patrik  Liba on 24/10/2025.
//

#ifndef CPPFM_BLACKSCHOLESFORMULAS_H
#define CPPFM_BLACKSCHOLESFORMULAS_H

#include "FinancialInstrument.h"
#include "DiscountCurve.h"
#include <cmath>

/**
 * @class BlackScholesFormulas
 * @brief Static utility class for Black-Scholes formulas and Greeks
 */
class BlackScholesFormulas
{
public:
    // d1/d2 helpers
    static double d1(double spot, double strike, double rate, double volatility, double maturity);
    static double d2(double spot, double strike, double rate, double volatility, double maturity);
    static double d1(double spot, double strike, const DiscountCurve& discountCurve, double volatility, double maturity);
    static double d2(double spot, double strike, const DiscountCurve& discountCurve, double volatility, double maturity);
    
    // Pricing functions
    static double price(double spot, double strike, double rate, double volatility, double maturity, Option::Type optionType);
    static double callPrice(double spot, double strike, double rate, double volatility, double maturity);
    static double putPrice(double spot, double strike, double rate, double volatility, double maturity);
    static double price(double spot, double strike, const DiscountCurve& discountCurve, double volatility, double maturity, Option::Type optionType);
    static double callPrice(double spot, double strike, const DiscountCurve& discountCurve, double volatility, double maturity);
    static double putPrice(double spot, double strike, const DiscountCurve& discountCurve, double volatility, double maturity);
    
    // Greeks - First order
    static double delta(double spot, double strike, double rate, double volatility, double maturity, Option::Type optionType);
    static double delta(double spot, double strike, const DiscountCurve& discountCurve, double volatility, double maturity, Option::Type optionType);
    static double vega(double spot, double strike, double rate, double volatility, double maturity);
    static double vega(double spot, double strike, const DiscountCurve& discountCurve, double volatility, double maturity);
    static double theta(double spot, double strike, double rate, double volatility, double maturity, Option::Type optionType);
    static double theta(double spot, double strike, const DiscountCurve& discountCurve, double volatility, double maturity, Option::Type optionType);
    static double rho(double spot, double strike, double rate, double volatility, double maturity, Option::Type optionType);
    
    // Greeks - Second order
    static double gamma(double spot, double strike, double rate, double volatility, double maturity);
    static double gamma(double spot, double strike, const DiscountCurve& discountCurve, double volatility, double maturity);
    static double vanna(double spot, double strike, double rate, double volatility, double maturity);
    static double volga(double spot, double strike, double rate, double volatility, double maturity);
    
private:
    BlackScholesFormulas() = delete;
};


#endif //CPPFM_BLACKSCHOLESFORMULAS_H