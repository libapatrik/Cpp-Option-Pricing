//
// Created by Patrik  Liba on 17/09/2025.
//

#ifndef CPPFM_FINANCIALPAYOFFS_H
#define CPPFM_FINANCIALPAYOFFS_H

class FinancialInstrument // Abstract class - base class
{
public:

    virtual FinancialInstrument* clone() const = 0; // Virtual Pure method - clone

    virtual ~FinancialInstrument() = default;       // Virtual destructor for FinancialInstrument

};

class Option : public FinancialInstrument // Abstract class - base class for Option
{
public:
    enum class Type {Call, Put}; // place first!

    explicit Option(Option::Type type); // Constructor with parameter
    Option() = default;                 // Default constructor
    Option(const Option&) = default;
    Option& operator=(const Option&   ) = default;
    ~Option() override = default;       // ALWAYS DECLARE BASE CLASS DESTRUCTOR VIRTUAL -> Avoid memory leak

    // Getter as virtual - allow model to check if its call or put
    // TODO: Each Model or Option will override this?

    // inline virtual Type getType() const // was confusing, so used getType instead of Type()
    // {
    //     return _type; // pure virtual method
    // }
    virtual Type getType() const = 0; // pure virtual method - getter
    Option* clone() const override = 0; // Virtual Pure method - clone

protected:
    Type _type;

};


class EuropeanOptionPayoff : public Option
{
public:
    EuropeanOptionPayoff(Type type);               // Constructor with parameter - sets to either put or call when creating object
    EuropeanOptionPayoff(const EuropeanOptionPayoff&) = default;
    EuropeanOptionPayoff& operator=(const EuropeanOptionPayoff&) = default;
    ~EuropeanOptionPayoff() override = default;

    Option::Type getType() const override; // C/P provided by Option class; EUOption overrides and provides its own implementation
    EuropeanOptionPayoff* clone() const override;

private:
    // Type _type;
};




#endif //CPPFM_FINANCIALPAYOFFS_H
