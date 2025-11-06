/**
 * Learning the basics
 */
#include <gtest/gtest.h>  
#include <iostream>
/// PRACTICE: Simple excercise aside
struct BankAccount {
    double balance = 0;

    BankAccount() : balance(0.0) {}

    explicit BankAccount(const double balance)
        : balance(balance)
    {
    }

    void deposit(double amount)
    {
        balance += amount;
    }

    bool withdraw(double amount)
    {
        if (amount <= balance) {
            balance -= amount;
            return true;
        }
        return false;
    }
};

struct BankAccountTest : testing::Test
{
    BankAccount* account;

    BankAccountTest()
    {
        account = new BankAccount();
    }

    ~BankAccountTest()
    {
        delete account;
    }
};

struct accountState
{
    double initialBalance;
    double withdrawAmount;
    double finalBalance;
    bool success;

    // Friend function for printing accountState objects
    friend std::ostream& operator<<(std::ostream& os, const accountState& obj) {
        return os << "initial_balance: " << obj.initialBalance
                  << ", withdraw_amount: " << obj.withdrawAmount
                  << ", final_balance: " << obj.finalBalance
                  << ", success: " << obj.success;
    }
};

struct WithdrawAccountTest : BankAccountTest, testing::WithParamInterface<accountState>
{
    WithdrawAccountTest()
    {
        account->balance = GetParam().initialBalance;
    }
};

TEST_P(WithdrawAccountTest, FinalBalance)
{
    auto as = GetParam();
    auto success = account->withdraw(as.withdrawAmount);
    EXPECT_EQ(as.finalBalance, account->balance) << "Test parameters: " << as;
    EXPECT_EQ(as.success, success) << "Test parameters: " << as;
}

TEST_F(BankAccountTest, BankAccountStatsEmpty)
{
    EXPECT_EQ(0.0, account->balance); // New account should start with 0 balance
}

TEST_F(BankAccountTest, CanDepositMoney)
{
    account->deposit(100.0);
    EXPECT_EQ(100.0, account->balance);
}

INSTANTIATE_TEST_SUITE_P(Default, WithdrawAccountTest,
    testing::Values(accountState{100, 50, 50, true},
                    accountState{100, 200, 100, false}
                    ));

int main(int argc, char* argv[]) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}




// TODO: Include test_utils.h when ready

/*
PURPOSE: Test end-to-end workflows and realistic market scenarios
- Complete workflow: construct -> interpolate -> price
- Realistic market data scenarios
- Consistency between different methods
- Performance with large surfaces
- Builder pattern integration

TODO:
1. Test complete workflow
2. Test with realistic market data
3. Test consistency between methods
4. Test performance with large surfaces
5. Test builder integration
*/

// TODO: Add test cases