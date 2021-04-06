## Portfolio_Tracker

# Prerequisites necessary to use this code

# (1) Create an AlphaVantage account for free at  https://www.alphavantage.co/support/#api-key and request an API
#
# (2) Construct a spreadsheet (.csv) with  all sales/purchase information in it:
#
# The structure of the spreadsheet should be like this:
#
# Date       | Quantity_Stock1 | Price_Stock1    | Payments_Stock1| Cost_Stock1 |
# 20xx-xx-xx |      Float      |      Float      |      Float     |      Float  |
#
# "Stock1" should be the string-handle used to pull the security-series from AlphaVantage.
# E.g. if you want to pull Apple's stock price, use "Quantity_AAPL", "Price_AAPL", etc
#
# For symbol search define SEARCHWORD and find ticker:
# https://www.alphavantage.co/query?function=SYMBOL_SEARCH&keywords=SEARCHWORD&apikey=YOURKEY
#
# Adding new infomation just requires adding  additional rows. Adding a new security
# means adding four new columns
