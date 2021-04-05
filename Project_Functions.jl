function get_trading_info(wdir)

    # Function that loads all trading information available:
    #
    # Input: 'wdir' provides information on file_directory to load .csv
    #
    # Output:  - 'raw_trades' returns a dataframe of the .csv-file, allows to read out the data
    #          - 'trading_start' returns the date string of the first trade undertaken
    #          - 'header_types' returns header elements like "quantity", "price", "payments" etc.
    #          - 'name_list' returns names/ticker_symbols of all securities that are stored in file

    # Load trading information from file, needs to be provided by user.
    raw_trades  = DataFrame(CSV.File(wdir*"/trades.csv", datarow = 2))

    # Fill empty spaces/missing values with zeros
    raw_trades  = coalesce.(raw_trades, 0.0)

    # Find/sort headers

    headers       = string.(names(raw_trades)); # DataFrame header row
    header_types  = ["Quantity_", "Price_", "Payments_", "Cost_in_Euro_", "FX_USD_"]; # Different types of information, provided for each security
    n_titles      = (size(raw_trades,2) - 1)/length(header_types); # Determine the number of securities in spreadsheet

    # Every size(header_types,1)-step, there is information on a new security. Define index with spacing of interval
    var_pos       = convert.(Int64,hcat([(i-1)*length(header_types)+1 for i = 1:n_titles].+1)); # (+1) takes date-field on first position into account

    function get_ETF_names(var_pos, headers)

        # Function that, given the header row of the DataFrame, returns the names/ticker
        # symbols of securities. Information is needed to pull the market data.

        # Input:  - 'var_pos' is an index in which interval a new security arrives
        #         - 'headers' is a list with the header row of the dataframe
        # Output: - 'name_list' contains names/ticker symbols of securities.

        # predefine empty array
        name_list = Array{String,1}(undef,length(var_pos));

        # loop over index positions to extract names
        for (i,j) in enumerate(var_pos)
            length_ident = length("Quantity_")
            name_list[i] = headers[j][length_ident+1:end];
        end
        return name_list
    end

    # Get names of securities: Use this information to pull market data
    name_list     = get_ETF_names(var_pos, headers);
    # Data of first trade
    trading_start = string(raw_trades[1,:Date]);

    return raw_trades, trading_start, header_types, name_list
end

function get_initial_quantities(raw_trades,start_date, header_types)
    # Function that given trading information finds the initial quantities of
    # a security held in case 'start_date' is later than the 'trading_start'.
    # Input:  - 'raw_trades' A dataframe of the trading information.
    #         - 'header_types', all different information on a certain security.
    #         - 'start_date' is a Date string of the sample start.
    # Output: - 'initial_quantities' is an array of security quantities at trading start.
    #         - 'traded_quantities' is an binary array: either a security is traded within the sample or not.


    # find all titles in trading information
    headers       = string.(names(raw_trades));
    n_quantities  = convert(Int64,(size(raw_trades,2) - 1)/length(header_types));

    # add one for date-field on first position
    quantity_pos           = convert.(Int64,hcat([(i-1)*length(header_types)+1 for i = 1:n_quantities].+1));
    # find position until when you need to sum up traded quantities
    prev_trades_pos    = findall(x -> x >= Date(start_date), raw_trades.Date)[1]-1;
    trade_interval_pos = prev_trades_pos+1:size(raw_trades,1);

    function get_quantities(quantity_pos, headers,prev_trades_pos,trade_interval_pos)
        initial_quantities = Array{Float64,1}(undef,length(quantity_pos));
        traded_quantities   = Array{Float64,1}(undef,length(quantity_pos))
        for (i,j) in enumerate(quantity_pos)
            initial_quantities[i] = sum(raw_trades[1:prev_trades_pos, j]);
            traded_quantities[i]  = sum(abs.(raw_trades[trade_interval_pos[1]:trade_interval_pos[end], j]).> 0);
        end

        initial_quantities = abs.(round.(initial_quantities, digits=4))
        traded_quantities  = traded_quantities .> 0;

        return initial_quantities, traded_quantities
    end

    initial_quantities, traded_quantities   = get_quantities(quantity_pos, headers,prev_trades_pos,trade_interval_pos);

    return initial_quantities, traded_quantities
end

function rearrange_trade_data(traded_pos_start,traded_pos_end, header_types, name_list)

    # Function that constructs an index of securities actually traded in interval from raw_trades
    #
    # Input:  - 'traded_pos_start', 'traded_pos_end' are two indeces that mark start/end of all information on a certain security
    #         - 'header_types', all different information on a certain security.
    #         - 'name_list', (reduced) list of securities traded in selected horizon.
    # Output: - 'index_to_rearrange_raw_trades', index used to slice raw_trades into smaller dataframe.

    index_to_rearrange_raw_trades =  Array{Float64,1}(undef,length(header_types)*length(name_list)+1);
    index_to_rearrange_raw_trades[1] = 1;


    for (i,j) in enumerate(traded_pos_start)
        index_to_rearrange_raw_trades[(i-1)*length(header_types)+2:(((i+1)-1)*length(header_types)+1)] .= traded_pos_start[i]:traded_pos_end[i];
    end
    index_to_rearrange_raw_trades = convert.(Int64,index_to_rearrange_raw_trades)
    return  index_to_rearrange_raw_trades
end

function dynamic_tuple(name_list)
# Remove all information behind ticker (e.g. like .FRK or .DE) if there is any.
# Helps to use string more efficient in meta-programming.

# Input: 'name_list': contains trading /ticker symbol code
# Output: 'label':   size(name_list,1)-dimensional tuple of names attributed to traded securities.
# Example: 'AAPL' remains 'AAPL' but 'VOW3.DEX' would become only 'VOW3'.


    label = []
    for (i,j) in enumerate(name_list)
        j_label = Symbol(split(j,r"[.]")[1])
        label = [label..., j_label]
    end
    return label
end

function raw_to_dataframe(rawData)

    df = DataFrame(rawData[1])
    dfNames = Symbol.(vcat(rawData[2]...))
    df = DataFrames.rename!(df, dfNames)

    df.Date = Date.(df.timestamp)
    for x in (:open, :high, :low, :close, :adjusted_close, :dividend_amount)
        df[!, x] = Float64.(df[!, x])
    end
    df.volume = Int64.(df.volume)

    # Get rid of some entries not used
    select!(df, Not([:volume, :dividend_amount, :split_coefficient]))

    # add empty column for quantities held
    pos = size(df,2)+1
    insertcols!(df, pos, :Quantity_held => 0.0)
    insertcols!(df, pos, :Current_value => 0.0)
    insertcols!(df,pos,  :Invested_Amount => 0.0)

    return df
end

function intra_to_dataframe(rawData)

    df = DataFrame(rawData[1])
    dfNames = Symbol.(vcat(rawData[2]...))
    df = DataFrames.rename!(df, dfNames)

    df.DateTime = DateTime.(df.timestamp, "yyyy-mm-dd HH:MM:SS")
    for x in (:open, :high, :low, :close)
        df[!, x] = Float64.(df[!, x])
    end
    df.volume = Int64.(df.volume)

     # add new column for quantities held
     pos = size(df,2)+1
     insertcols!(df, pos, :Quantity_held => 0.0)
     insertcols!(df, pos, :Current_value => 0.0)
     insertcols!(df, pos, :Invested_Amount => 0.0)

    return df
end


function draw_series(name::String)
# Function that dowloads data using a ticker string-input
# Input:   - 'name' is the ticker string that matches the series to be obtained
# Output:  - 'df' of price information on the security.

    println("draw $name")
    raw = AlphaVantage.time_series_daily_adjusted(name, outputsize="full", datatype="csv");
    # reaarrange csv-data into DataFrame
    df   = raw_to_dataframe(raw);
    # Pause: if multiple series are downloaded, sleep makes sure you don't make more than 5 requests per minute,
    # as availabe in free API Version.

    sleep(8)

    start_sample = findall(x -> x == start_date, df.timestamp)[1];

    if isempty(end_date)
        end_sample = 1
    else
        end_sample   = findall(x -> x == end_date, df.timestamp)[1];

    end

    df = df[end_sample:start_sample,:];


    # Set initial quantities held
    global init_name_quantities
    df[!,:Quantity_held].= init_name_quantities[name]
    df[!,:Current_value].= df.Quantity_held.*df.adjusted_close

    return df
end

# Functions used to trim the data to a common length

function find_and_select(min_date,df,dates,key)
    """
    min_date = minimal length of all stocks in sample
    df = dataframe to be trimmed
    dates = named tuple of date vectors
    key = the key at which the named tuple is accessed (e.g. :tech, :emm, ...)
    """
    f = findall(x -> x in min_date, dates[key])
    return df[key][[f,1][1],:]
end


check_size(s) = size(unique(s))[1]


function same_size_df(data)
    # Find smallest dataframe
    dates  = [string.(s.timestamp) for s in data]
    size_selector = [size(dates[i])[1] for i in 1:length(dates)]
    dates = NamedTuple{Tuple(labels)}(dates) #construct a named tuple of timestamp vectors

    min_pos       = findall(x -> x == minimum(size_selector), size_selector) #unchanged
    min_size      = minimum(size_selector)
    min_date      = dates[min_pos[1]]

    out = map(k -> find_and_select(min_date,data,dates,k), labels) #Map over all labels [:heath,...]
     #create a named tuple for the trimmed data
    new_data      = NamedTuple{Tuple(labels)}(out)
    #After trimming, check the size of the new vector
    dates  = [string.(s.timestamp) for s in data]
    size_selector  = [size(dates[h])[1] for h in 1:(length(dates))]
    cond = check_size(size_selector) == 1 #check equal sizing and store result
    return (cond = cond, data = new_data)
end

function trim_dataframe!(data::NamedTuple)

    """
    Trims the dataframe. Check the initial size and store initial result as init.
    Only requires data as input
    """
    dates = [string.(s.timestamp) for s in data]
    size_selector = [size(dates[i])[1] for i in 1:length(dates)]

    init = check_size(size_selector) == 1
    if init
        println("Dataframes don't need to be adjusted. They already have the same size.")
    else
        satisfied = false #iterate while condition is not satisifed
        counter = 1
        while !satisfied
            println("Trimming $counter")
            adjusted = same_size_df(data)
            satisfied = adjusted.cond #assign result of size check to "satisfied"
            data = adjusted.data
            counter += 1
        end
    end
    return data
end

function recover_missing_data!(untrimmed, data, missing_days, trading_days,labels)

# Auxillary stuff
df_list  = string.(labels);
var_list      = [df_list[i]*"_ETF" for i = 1:length(labels)];
timeline = eval(Meta.parse("data."*string(df_list[1])*"[!,:Date]"))


    for i = 1:length(df_list)

        trimmed_aux      = eval(Meta.parse("data."*string(df_list[i])))
        untrimmed_aux    = eval(Meta.parse("untrimmed."*String(df_list[i])))
        missing_position = findall(x -> x in missing_days, untrimmed_aux[!,:Date])

            if isempty(missing_position)
                missing_days_update = missing_days
                while isempty(missing_position)
                    missing_days_update  = missing_days_update-Dates.Day(1)
                    missing_position     = findall(x -> x in missing_days_update, untrimmed_aux[!,:Date])
                end
                # set to value of last available trading day
                missing_data     = untrimmed_aux[missing_position,:]
                # overwrite date values --> needed for merge to work
                missing_data[!,:Date]      = missing_days
                trimmed_update = join(missing_data, trimmed_aux, on = names(trimmed_aux), kind = :outer)
                trimmed_update   = sort(trimmed_update, (:Date), rev = true)

                data = merge(data, (keys(data)[i] => trimmed_update,))

            elseif length(missing_position) < length(missing_days)

                missing_data     = untrimmed_aux[missing_position,:]
                trimmed_update   = join(missing_data, trimmed_aux, on = names(trimmed_aux), kind = :outer)
                trimmed_update   = sort(trimmed_update, (:Date), rev = true)

                remaining_days   = findall(x -> x == 1, missing_days .∉ [trimmed_update[!,:Date]])
                remained_missing = missing_days[remaining_days]

                remaining_position = findall(x -> x in remained_missing, untrimmed_aux[!,:Date])

                remaining_days_update = remained_missing
                while isempty(remaining_position)
                    remaining_days_update   = remaining_days_update - Dates.Day(1)
                    remaining_position      = findall(x -> x in remaining_days_update, untrimmed_aux[!,:Date])
                end

                # set to value of last available trading day
                missing_data     = untrimmed_aux[remaining_position,:]
                # overwrite date values --> needed for merge to work
                missing_data[!,:Date]      = remained_missing

                trimmed_update = join(missing_data, trimmed_update, on = names(trimmed_update), kind = :outer)
                trimmed_update   = sort(trimmed_update, (:Date), rev = true)

                data = merge(data, (keys(data)[i] => trimmed_update,))

            else
                missing_data     = untrimmed_aux[missing_position,:]
                trimmed_update   = join(missing_data, trimmed_aux, on = names(trimmed_aux), kind = :outer)
                trimmed_update   = sort(trimmed_update, (:Date), rev = true)

                data = merge(data, (keys(data)[i] => trimmed_update,))

            end
    end

    return data
end

function calc_portfolio_stats(var_list, df_list)

# Total Portfolio Value
        for (i,j) in enumerate(var_list)
            rhs_value  = var_list[i]*"[!,:Current_value]";
            if i == 1
                lhs_value = "total_value =";
            else
                lhs_value = "total_value +="
            end
            iter_string_value = Meta.parse(lhs_value*rhs_value)
            eval(iter_string_value)
        end

# Total Invested_Amount
        for (i,j) in enumerate(var_list)
            rhs_value  = var_list[i]*"[!,:Invested_Amount]";
            if i == 1
                lhs_value = "total_amount_invested =";
            else
                lhs_value = "total_amount_invested +="
            end
            iter_string_value = Meta.parse(lhs_value*rhs_value)
            eval(iter_string_value)
        end


# Portfolio Weights

        for (i,j) in enumerate(var_list)
            rhs_weights  = var_list[i]*"[!,:Current_value]./total_value*100";
            if i == 1
                lhs_weights = "portfolio_weights =";
                iter_string_weights = Meta.parse(lhs_weights*rhs_weights)
                eval(iter_string_weights)
            else
                lhs_weights = "weights =";
                iter_string_weights = Meta.parse(lhs_weights*rhs_weights)
                eval(iter_string_weights)
                hcat_string = Meta.parse("portfolio_weights = hcat(portfolio_weights, weights)")
                eval(hcat_string)
            end
        end

    # Calculate total return
    total_return = (total_value./total_amount_invested.-1).*100;


    return total_value, total_amount_invested, total_return, portfolio_weights
end


function Bollinger_Bands(data,T_Roll,MA_window)

    mean_var = zeros(T_Roll,1);
    std_var  = zeros(T_Roll,1);

    for i_Boll = 1:T_Roll
        mean_var[i_Boll] = mean(data[i_Boll:MA_window+i_Boll-1,:adjusted_close]);
        std_var[i_Boll] = std(data[i_Boll:MA_window+i_Boll-1,:adjusted_close]);
    end

    Bollinger_Bands = hcat(mean_var.-2*std_var,data[MA_window+1:end,:adjusted_close],mean_var, mean_var.+2*std_var)
    return Bollinger_Bands
end
