% Importing_lecture

% load a previously saved .mat file


% Excel importl
excel_import = xlsread('excel_random_numbers.xlsx') 

% or import a subset

excel_import = xlsread('excel_random_numbers.xlsx', 'a1:b2')

% CSV aka Comma Separated Values

csv_data = csvread('csv_random_numbers.csv')

% to import a subset ofo CSV, need to know that the top left value in the
% file is (0,0)

csv_data_subset = csvread('csv_random_numbers.csv', 0,0, [0 0 1 1])
% filename, first row, first column, [first row first col last row last col


% more stuff on ofunctions

% transpose; flip around an array

new_csv = transpose(csv_data_subset)

doc transpose

%size tells you the size of an array
% helps you avoid hard-coding numbers into your code. 

% bad don't do this

new_array = zeros(5,10); 


% better way:
num_of_rows = 5;
num_of_columns = 10'
new_array = zeros(num_of_rows, num_of_columns);
second_array = zeros(num_of_rows, num_of_coluumns);


% even better way:
new_array = zeros(num_of_rows, num_of_columns);
secondarray_ = zeros(size(new_array)); 

%statistics

sv_value = mean(csv_data_subset)
% mean gives mean of each column

% if you want the average of all elements;
 % : means "all" 

av_value = mean(csv_data_subset(:)) 


% Max value of an array:
highest_value = max(csv_data_subset(:))

% Max value within column 1 only:
highest_value_first_column = max(csv_data_subset(:,1))

% Of the second column 
highest_value_second_column = max(csv_data_subset(:,2))

% Make an array with summary statistics

% let row 1 be the mean 
% row 2 = median
% row 3 = first quartile (25th percentile)
% row 4 = third quaurtile (75th percentile)

summary_csv = zeros(4,1);

summary_csv(1) = mean(csv_data_subset(:));
summary_csv(2) = prctile(csv_data_subset(:), 50);
summary_csv(3) = prctile(csv_data_subset(:), 25);
summary_csv(4) = prctile(csv_data_subset(:), 75);