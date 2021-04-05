cd("File Path")
wdir     = pwd()
# # Write API Key to environment

ENV["ALPHA_VANTAGE_API_KEY"] = "YOUR KEY HERE"

# Core Packages needed
using AlphaVantage
using DataFrames
using DataFramesMeta
using Dates
using CSV

# For making nice plots
using Plots
using LaTeXStrings
using Statistics
using ColorSchemes
using TimeSeries

include("Project_Functions.jl")

################################################################################
                                    # Main Part #
################################################################################

## Prerequisites necessary to use this code

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
# https://www.alphavantage.co/query?function=SYMBOL_SEARCH&keywords=SEARCHWORD&apikey=6M3D8EV1HA6WKQFI
#
# Adding new infomation just requires adding  additional rows. Adding a new security
# means adding four new columns

## Procedures

# Set start/end date of of data sample:

start_date          = "2019-01-04";  # Trading information available starts on "2018-04-03".
end_date            = [];  # if "end_date" is empty, the last available value will be set.


# (1) Load Trading Information avialable: - Which titles and in which quantity have been traded on this account over the entire horizon?

raw_trades, trading_start, header_types, full_name_list = get_trading_info(wdir);

# (2) Use available trading information and match it to desired sample horizon.
#     Only load security data of securities traded in respective interval.
#
#     This requires the following steps:

# (2a) - Which titles and which quantity did I trade within selected start/end_date (today)?
#      - Find securities in "full_name_list" that are not in trading sample within selected time span.
#      - Create a subset that only pulls, analyzes and prints securities available in in
#        the selected time span.
# (2b) - Determine initial/already exisiting positions of a security for any date
#        starting after 'trading_start'. Neded for accounting/valuation purposes later.

if Date(start_date) < Date(trading_start)

    # Date to start portfolio analysis: should not set to values before 'trading_start'
    start_date           = trading_start;
    println("Warning: 'start_date' set to "*trading_start*" , as quantities held are zero before.")

    # Initial quantities are zero by definition
    initial_quantities, traded_securities  =  get_initial_quantities(raw_trades,start_date, header_types);
    init_name_quantities = Dict(full_name_list[i] => initial_quantities[i] for i = 1:length(initial_quantities));

    # Find which funds where traded in regarded time period
    name_list = full_name_list[traded_securities.> 0];

    # Extract securities actually traded in interval from raw_trades: Save in sub-DataFrame 'trades'.
    # If full horizon is selected, all securities and all trades are selected. 'trade' will be identical to 'raw_trades'.

    traded_securities = findall(x -> x in name_list,full_name_list);
    traded_pos_start  = convert.(Int64,hcat([(j-1)*length(header_types)+1 for (i,j) in enumerate(traded_securities)].+1));
    traded_pos_end    = traded_pos_start.+(length(header_types)-1);

    index_to_rearrange_raw_trades = rearrange_trade_data(traded_pos_start,traded_pos_end, header_types, name_list);
    trades = raw_trades[:,index_to_rearrange_raw_trades]

elseif Date(start_date) >= Date(trading_start)

    # Find initial quantities, update raw_trades to have only those securites that were traded
    initial_quantities, traded_securities    =  get_initial_quantities(raw_trades,start_date, header_types);
    init_name_quantities = Dict(full_name_list[i] => initial_quantities[i] for i = 1:length(initial_quantities))

    # Find which funds where traded in regarded time period
    name_list = full_name_list[traded_securities.> 0];

    # Extract securities actually traded in interval from raw_trades: Save in sub-DataFrame 'trades'.

    traded_securities = findall(x -> x in name_list,full_name_list);
    traded_pos_start  = convert.(Int64,hcat([(j-1)*length(header_types)+1 for (i,j) in enumerate(traded_securities)].+1));
    traded_pos_end    = traded_pos_start.+(length(header_types)-1);

    index_to_rearrange_raw_trades = rearrange_trade_data(traded_pos_start,traded_pos_end, header_types, name_list);
    trades = raw_trades[:,index_to_rearrange_raw_trades]

else
    println("Warning: Something went wrong when defining the date value. Check the spreadsheet.")

end

# (3) Obtain price information and clean raw data

# (3a) Preparations

# Labels of reduced security list in symbol notation (used for meta-programming)
labels      = dynamic_tuple(name_list)

# Labels of full security list in symbol notation (used for meta-programming)
full_labels = dynamic_tuple(full_name_list);

# String version of labels: Used as identifier in loops and to access list of variables
df_list     = string.(labels);
var_list    = [df_list[i]*"_Stock" for i = 1:length(labels)]

# Produce header row in symbol notation
csv_headers = reduce(vcat,[header_types.*df_list[i] for i = 1:length(labels)]);
# Add first column header entry
csv_headers = ["Date"; csv_headers];

csv_headers_full = reduce(vcat,[header_types.*string.(full_labels)[i] for i = 1:length(full_labels)]);
# Add first column header entry
csv_headers_full = ["Date"; csv_headers_full];

# Exchange headers in dataframes such that they are accessible with meta-programming using symbol notation.
trades        = names!(trades, Symbol.(csv_headers));
raw_trades    = names!(raw_trades, Symbol.(csv_headers_full))

# (3b) Download Market Data

# Mapping that applies the "draw_series"-function to every variable in name_list

# Raw data: It can happen that dataframes are of different length (holidays).
# We will take care of that in following steps.
vals          = map(draw_series,name_list);

# Construct NamedTuples of raw data, which makes data moer accessible.
# i.e. untrimmed.'label[1]' returns a DataFrame of the first security.
untrimmed     = NamedTuple{Tuple(labels)}(vals);

# Get all stocks to the same size (Differences due to different trading days)

# Trim data to same size
data          = trim_dataframe!(untrimmed);

# Export securities as individual DataFrames

for (i,j) in enumerate(labels)
    rhs = "data["*string(i)*"]";
    lhs = String(j)*"_Stock";
    eval(Meta.parse(lhs*"="*rhs))
end


# (4) Adjust portfolio information for trades and check if any day at which trading took place
# has been deleted in trimming process.

# Define Timeline: Choose whether stock/ETF you like as all df should be of the same length
timeline      = data[1].Date

start_index   = minimum(findall(x -> x >= Date(start_date), raw_trades.Date));
if isempty(end_date)
 end_index = maximum(findall(x -> x >= Date(start_date), raw_trades.Date));
else
    end_index     = maximum(findall(x -> x <= Date(end_date), raw_trades.Date));
end

# # Further adjust dataframe with trading info accordingly
trades        = trades[start_index:end_index,:];

# Find all days in which sales / purchases took place
trading_days  = trades[!, :Date]
trade_index   = sort(findall(x -> x in trading_days, timeline), rev = true)
#
trading_days_in_timeline = timeline[trade_index]


##
if trading_days == trading_days_in_timeline

else

    println("Warning: There are days at which trades took place that have been thrown out of the sample. Think about altering the trimming routine")

    missing_days = findall(x -> x == 1, trading_days .∉ [trading_days_in_timeline])
    missing_days = trading_days[missing_days]

    data = recover_missing_data!(untrimmed, data, missing_days,trading_days,labels)

    # Export (again) as individual DataFrames
    for (i,j) in enumerate(labels)
        rhs = "data["*string(i)*"]";
        lhs = String(j)*"_Stock";
        eval(Meta.parse(lhs*"="*rhs))
    end

    # Reset Timeline
    timeline      = data[1].Date

    trade_index   = sort(findall(x -> x in trading_days, timeline), rev = true)

    trading_days_in_timeline = timeline[trade_index]

    ##
    if trading_days == trading_days_in_timeline
    println("Discarded information on trading day data has been succesfully recovered. Missing data has been set to the last available trading day.")
    else
    println("Something went wrong. Check the algorithm for potential undetected bugs.")
    end
end


##

for (i,v) in enumerate(df_list)

    # add quantities
    df_new_entry = Meta.parse("Cum_Quantity_"*v)
    df_entry     = Meta.parse("Quantity_"*v)
    #raw_trades[!,df_new_entry].= initial_quantities[i].+cumsum(raw_trades[!,df_entry])
    raw_trades[!,df_new_entry].= cumsum(raw_trades[!,df_entry])

    # add information about invested amount, payments and cost, price_tag
    df_invested_amount = Meta.parse("Invested_Amount_"*v);
    df_takeouts        = Meta.parse("Payouts_"*v);
    df_payments_entry  = Meta.parse("Payments_"*v);
    df_cost_entry      = Meta.parse("Cost_in_Euro_"*v);
    df_price_tag       = Meta.parse("Price_"*v);

    raw_trades[!,df_payments_entry] .= raw_trades[!,df_payments_entry];
    raw_trades[!,df_cost_entry]     .= raw_trades[!,df_cost_entry];
    raw_trades[!,df_price_tag]      .= raw_trades[!,df_price_tag];
    raw_trades[!,df_invested_amount].= 0.0;

    # Calculate the amount invested and adjust accordingly it after a sale occurs (Fifo method)
    for i_trade = 1:size(raw_trades,1)

    # special condition of first entry of Invested Amount (nothing to cumulate)

        if i_trade == 1
            raw_trades[i_trade, df_invested_amount] = raw_trades[i_trade,df_entry]*raw_trades[i_trade,df_price_tag];
        else

            if raw_trades[i_trade,df_entry] >= 0
                raw_trades[i_trade, df_invested_amount] = raw_trades[i_trade-1, df_invested_amount]+raw_trades[i_trade,df_entry]*raw_trades[i_trade,df_price_tag];
            elseif raw_trades[i_trade,df_entry] < 0

                amount_to_sell         = raw_trades[i_trade,df_entry];
                #cumulated_quantity_aux = raw_trades[!,df_new_entry];
                cumulated_quantity_aux = cumsum(raw_trades[!,df_entry]);
                remaining_quantity     =  cumulated_quantity_aux[i_trade-1] - abs(amount_to_sell); # cumulated quantity up to current sale

                if remaining_quantity > 0 # not everything has been sold -> calculate remaining invested amount according to FIFO accounting

                    # Find ALL purchases and constrain index to all purchases BEFORE current sale (and even time horizon considered)


                    all_purchases =  reverse(findall( x -> x > 0, raw_trades[!,df_entry]));
                    buy_index     =  all_purchases[all_purchases.<i_trade];

                    if isempty(buy_index)
                    all_purchases =  reverse(findall( x -> x > 0, raw_trades[!,df_entry]));
                    buy_index     =  all_purchases[all_purchases.<(i_trade)];

                    else
                    end

                    check_aux = remaining_quantity;
                    i_sell    = 0;

                    while check_aux > 0

                        i_sell = i_sell +1;

                        if raw_trades[buy_index[i_sell],df_entry] < check_aux #trades.Quantity_Tech(buy_index(i_sell)) < check_aux

                            # Case: Last purchase was smaller than overall
                            # quantity remaining. Remaining amount has to be
                            # allotted the the previous purchases until every
                            # part of the remaining amount has been allocated.

                            raw_trades[i_trade, df_invested_amount] = raw_trades[i_trade, df_invested_amount] + raw_trades[buy_index[i_sell], df_entry]*raw_trades[buy_index[i_sell], df_price_tag];
                            check_aux = cumulated_quantity_aux[buy_index[i_sell+1]] - abs(raw_trades[i_trade, df_entry]);

                        elseif raw_trades[buy_index[i_sell],df_entry] > check_aux

                            # Case: Last purchase was higher than overall
                            # quantity remaining. Invested amount has to be
                            # priced using the price of the last purchase and
                            # the remaining quantity.
                            raw_trades[i_trade, df_invested_amount] = raw_trades[i_trade, df_invested_amount] + check_aux*raw_trades[buy_index[i_sell], df_price_tag];
                            check_aux = 0; # used to exit while-condition - has to hold in this case as all remaining quantity is priced with price of the last purchase

                        end

                    end

                else # everything has been sold -> remaining invested amount has to be zero
                end
            else
            end
        end
    end




end

##

trades        = raw_trades[start_index:end_index,:];

for (i,v) in enumerate(df_list)

    # add quantities
    df_new_entry = Meta.parse("Cum_Quantity_"*v)

    # add information about invested amount, payments and cost, price_tag
    df_invested_amount = Meta.parse("Invested_Amount_"*v);

    # Broadcast Variables to full dataframe over all dates
    for (j,k) in enumerate(trade_index)
        # Quantity
        lhs = var_list[i]*"["*string(1:k) *",:Quantity_held]"
        rhs = " .= trades["*string(j)*",:"*string(df_new_entry)*"]"
        iter_string = Meta.parse(lhs*rhs)
        eval(iter_string)

        # Invested_Amount
        lhs_invs_ = var_list[i]*"["*string(1:k) *",:Invested_Amount]"
        rhs_invs_ = " .= trades["*string(j)*",:"*string(df_invested_amount)*"]"
        iter_string_invs = Meta.parse(lhs_invs_*rhs_invs_)
        eval(iter_string_invs)
    end

    lhs = var_list[i]*"[!,:Current_value]"
    rhs = " .="*var_list[i]*"[!,:adjusted_close].*"*var_list[i]*"[!,:Quantity_held]"
    iter_string    = Meta.parse(lhs*rhs)
    eval(iter_string)

end



# (5) Calculate some portfolio stats: Total Value, Returns, Portfolio Weights etc
total_value, total_amount_invested, total_return, portfolio_weights = calc_portfolio_stats(var_list, df_list)
