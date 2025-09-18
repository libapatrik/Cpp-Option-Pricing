//
// Created by Patrik  Liba on 17/09/2025.
//

#ifndef CPPFM_FINANCIALPAYOFFS_H
#define CPPFM_FINANCIALPAYOFFS_H

class FinancialPayoff
{
    virtual ~FinancialPayoff() = default; // Virtual destructor
};

class Option : public FinancialPayoff
{
public:
    // Default constructor
    Option() = delete;
    // Constructor with parameters
    Option(double spot, double mu, double sigma);
    // Copy constructor
    Option(const Option& model);
    // Clone method
    Option* clone() const override;

    // Copy Assignment Operator
    Option& operator=(const Option& model);

    // Destructor
    ~Option() override = default; // ALWAYS DECLARE BASE CLASS DESTRUCTOR VIRTUAL -> Avoid memory leak


    enum class Type { Call, Put };

    // Getter as virtual - allow model to check if its call or put
    //TODO: Each Model or Option will override this?
    inline virtual Type getType() const // was confusing, so used getType instead of Type()
    {
        return _type; // pure virtual method
    }

protected:
    Type _type;

};

class EuropeanOptionPayoff : public Option
{
public:
    EuropeanOptionPayoff(Type type);


    Type getType() const override; // C/P provided by Option class; EUOption overrides and provides its own implementation


private:
    // Type _type;
};


#endif //CPPFM_FINANCIALPAYOFFS_H