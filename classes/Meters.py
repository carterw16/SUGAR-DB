"""The class for meters for Three-phase power flow
File Created: 04/12/2017"""


class Meters:

    def __init__(self,
                 name,
                 ID,
                 phases,
                 nominal_voltage,
                 measured_real_energy,
                 measured_reactive_energy,
                 measured_power,
                 measured_powerA,
                 measured_powerB,
                 measured_powerC,
                 measured_demand,
                 measured_real_power,
                 measured_reactive_power,
                 measured_voltageA,
                 measured_voltageB,
                 measured_voltageC,
                 measured_currentA,
                 measured_currentB,
                 measured_currentC,
                 bill_day,
                 price,
                 monthly_fee,
                 monthly_bill,
                 previous_monthly_bill,
                 monthly_energy,
                 previous_monthly_energy,
                 power_market,
                 first_tier_price,
                 second_tier_price,
                 third_tier_price,
                 first_tier_energy,
                 second_tier_energy,
                 third_tier_energy,
                 bustype=2,
                 bill_mode=None):
        self.name = name
        self.ID = int(ID)
        self.phases = int(phases)
        self.nominal_voltage = int(nominal_voltage)
        self.bustype = int(bustype)
        #
        self.measured_real_energy = measured_real_energy
        self.measured_reactive_energy = measured_reactive_energy
        self.measured_power = measured_power
        self.measured_powerA = measured_powerA
        self.measured_powerB = measured_powerB
        self.measured_powerC = measured_powerC
        self.measured_demand = measured_demand
        self.measured_real_power = measured_real_power
        self.measured_reactive_power = measured_reactive_power
        self.measured_voltageA = measured_voltageA
        self.measured_voltageB = measured_voltageB
        self.measured_voltageC = measured_voltageC
        self.measured_currentA = measured_currentA
        self.measured_currentB = measured_currentB
        self.measured_currentC = measured_currentC
        self.bill_day = bill_day
        self.price = price
        self.monthly_fee = monthly_fee
        self.monthly_bill = monthly_bill
        self.previous_monthly_bill = previous_monthly_bill
        self.previous_monthly_energy = previous_monthly_energy
        self.bill_mode = bill_mode
        self.monthly_energy = monthly_energy
        self.power_market = power_market
        self.first_tier_price = first_tier_price
        self.second_tier_price = second_tier_price
        self.third_tier_price = third_tier_price
        self.first_tier_energy = first_tier_energy
        self.second_tier_energy = second_tier_energy
        self.third_tier_energy = third_tier_energy
