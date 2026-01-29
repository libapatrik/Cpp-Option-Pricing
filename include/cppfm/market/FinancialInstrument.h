//
// Created by Patrik  Liba on 17/09/2025.
//

#ifndef CPPFM_FINANCIALPAYOFFS_H
#define CPPFM_FINANCIALPAYOFFS_H

/* TODO:
    Options:
        American Option -> Pricing via finite difference (Crank-Nicholson)
            Can try implement PDEs in general -> implicit, explicit, CN -> do generic formula and then per method add the correct weights
        Asian Option -> path dependent
    Other asset classes:
        Variance swaps
*/


class FinancialInstrument // Abstract class - base class
{
public:

    // payoff method -> value that you get at "expiry" -> something to think of

    virtual FinancialInstrument* clone() const = 0; // Virtual Pure method - clone

    virtual ~FinancialInstrument() = default;       // Virtual destructor for FinancialInstrument

};

class Option : public FinancialInstrument // Abstract class - base class for Option
{
public:
    enum class Type { Call, Put }; // place first!

    explicit Option(Option::Type type); // Constructor with parameter
    Option() = default;                 // Default constructor
    Option(const Option&) = default;
    Option& operator=(const Option&   ) = default;
    ~Option() override = default;       // ALWAYS DECLARE BASE CLASS DESTRUCTOR VIRTUAL -> Avoid memory leak

    // Getter as virtual - allow model to check if its call or put
    // Getter - accessor 
    virtual Type type() const = 0; // pure virtual method - simple accessor
    Option* clone() const override = 0; // Virtual Pure method - clone

protected:
    Type _type;
    
};


class EuropeanOptionPayoff : public Option
{
public:
    EuropeanOptionPayoff(Type type);               // Constructor with parameter - sets to either put or call when creating object
    EuropeanOptionPayoff(Type type, double strike, double maturity);  // Constructor with parameter - sets to either put or call when creating object
    EuropeanOptionPayoff(const EuropeanOptionPayoff&) = default;
    EuropeanOptionPayoff& operator=(const EuropeanOptionPayoff&) = default;
    ~EuropeanOptionPayoff() override = default;

    Option::Type type() const override; // C/P provided by Option class; EUOption overrides and provides its own implementation
    EuropeanOptionPayoff* clone() const override;

    // Virtual getters for polymorphic access 
    virtual double strike() const { return _strike; }
    virtual double maturity() const { return _maturity; }

private:
    // Option is specified by strike and maturity
    double _strike {0.0};
    double _maturity {0.0};
};


class AmericanOption : public Option
{
public:

};


#endif //CPPFM_FINANCIALPAYOFFS_H
