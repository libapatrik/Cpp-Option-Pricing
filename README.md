## UML Class Diagram

![UML Diagram](//www.plantuml.com/plantuml/png/XLN1Rjim3BthAtJRoN-m3LtR5Hc680Do77PGP2Oc8ai6HQ_8q8yVN2jkn4YmJz7Z8wd7HySFaGtns3l5OT0Sz24mXrwtI60PWrVnZ7umP7hjq0DV1j21k0jAVT9VfaTlGmjERn38aEqkCwhC8J8Pozew4-DCKlFlHsViUM-xu23Uui_ZW0HpEPsm3KGEVcQfkqJTD6zVK-qiUtfkWwbCCJWfUGbaGVy5OhJ6wSZ6dj6Vz35GqdUYC0ugTVtTOTk2_zE6LzyKz4X4fYIfykAIbgbVvls2tdYBnkaVFZADFNXU6_-Iyv4YARlXtU_vy7saC83n6VYj1q8hmOW2sN47zD4dAbNehZ-_DSsrPpqkI-NtcarCHRBRYL341TSw2wruNxXPseHIHBYtg5jdIfrpggG9L_AZ90ktygCWPoukY_QBq4w3s09teLrRiSIUoUhwrdk69qCkfAXQrX4Pdl5u109l8Oy-xVPpbIkixkhSlCIYSQLqYIhxw7hLbJN9o6b1CapCkp2OEqIcr6PMmMMw57jZi-5Exx4VdU50QpPFx8My_5iueAsXyM9wKo9iqKK6zANO8mNJ8lmq5v6gwcMJHHczs6CnLB9hNSli5VipSh8_nuGPTMLrqR98PQX5LAPMPxuinIOG_WkU9wc3UL-e0_0hWDTTMX24E9FbLv7mr75J_EOKC9adNYQjyhDWIwCUmBTZx_u3)

```plantuml
@startuml
top to bottom direction
skinparam linetype ortho

' ===== Model Hierarchy =====
abstract class Model
class BlackScholesModel
class DupireModel  
class HestonModel

BlackScholesModel --|> Model
DupireModel       --|> Model
HestonModel       --|> Model

' Force vertical stacking for Model group
Model -[hidden]down- BlackScholesModel
BlackScholesModel -[hidden]down- DupireModel
DupireModel -[hidden]down- HestonModel

' ===== Financial Instrument Hierarchy =====
abstract class FinancialInstrument
abstract class Option
class EuropeanOption
class AmericanOption
class AsianOption
class OtherInstruments

FinancialInstrument <|-- Option
FinancialInstrument <|-- OtherInstruments
EuropeanOption --|> Option
AmericanOption --|> Option
AsianOption    --|> Option

' Force vertical stacking for FinancialInstrument group
FinancialInstrument -[hidden]down- Option
FinancialInstrument -[hidden]down- OtherInstruments
Option -[hidden]down- EuropeanOption
EuropeanOption -[hidden]down- AmericanOption
AmericanOption -[hidden]down- AsianOption

' ===== Pricing Method Hierarchy =====
abstract class PricingMethod
class MonteCarloPricing
class FiniteDifferencePricing
class COSMethodPricing

MonteCarloPricing       --|> PricingMethod
FiniteDifferencePricing --|> PricingMethod
COSMethodPricing        --|> PricingMethod

' Force vertical stacking for PricingMethod group
PricingMethod -[hidden]down- MonteCarloPricing
MonteCarloPricing -[hidden]down- FiniteDifferencePricing
FiniteDifferencePricing -[hidden]down- COSMethodPricing

' ===== Path Simulator Hierarchy =====
abstract class PathSimulator
class EulerMaruyamaPathSimulator
class MilsteinPathSimulator
class BroadieKayaPathSimulator

EulerMaruyamaPathSimulator --|> PathSimulator
MilsteinPathSimulator      --|> PathSimulator
BroadieKayaPathSimulator   --|> PathSimulator

' Force vertical stacking for PathSimulator group
PathSimulator -[hidden]down- EulerMaruyamaPathSimulator
EulerMaruyamaPathSimulator -[hidden]down- MilsteinPathSimulator
MilsteinPathSimulator -[hidden]down- BroadieKayaPathSimulator

' ===== Horizontal spacing between hierarchies =====
Model -[hidden]right- FinancialInstrument
FinancialInstrument -[hidden]right- PricingMethod  
PricingMethod -[hidden]right- PathSimulator

@enduml
```
