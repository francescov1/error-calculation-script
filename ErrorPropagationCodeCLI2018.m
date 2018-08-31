%% Uncertainty Calculation
% Created by Francesco Virga, May 2017

% Simplified process:
% - get function from user
% - extract variable names from function
% - get values and errors for each variable from user
% - computer final values and partial derivatives
% - compute uncertainty contributions and plot

% Code uses predefined error messages when errors occur to simplify debugging. 
% If more in depth debugging is required, refer to the 'exception' variable in workspace

% Useful commands: 
% - 'clear': remove all variables in workspace 
% - 'clc': clear command window (does not clear workspace)
% - 'CTRL+C': terminates running program 

input(['\nUncertainty calculations using formula ''yErr^2 = (dy/dx1 * x1Err)^2 + (dy/dx2 * x2Err)^2 + ...''\n' ... % initial prompt
    '***If error occurs and more debugging is required, refer to exception variable in workspace***\n'...
    'Press enter to begin\n'],'s');


%% Get function
% - This section gets the function from the user that will be used for the calculation 

if exist('f','var') && exist('functionString','var') % if function already exists in workspace
    yn = input(sprintf('To use function ''%s'' again press Enter, otherwise type ''n''\n', functionString),'s');
    if strcmp(yn,'n') % if user wants to specify new equation
        [f, functionString, desiredVariable, variableNames, n] = getFunctionFromUser();
    end
else % if function does not exist in workspace
   [f, functionString, desiredVariable, variableNames, n] = getFunctionFromUser();
end


%% Get data
% - This section gets input variable values and their respective errors
% from the user

% initialize arrays
x = zeros(1, n);
dx = zeros(1, n);

% loop through variable names, prompting user to enter their values and
% errors
for i=1:n
    x(i) = input(sprintf('%s = ', variableNames(i)));
    dx(i) = input(sprintf('error in %s = ', variableNames(i)));
end

%% Compute value +/- uncertainty and plot
% - This section calculates the final value of the desired variable, its uncertainties and each of their contributions.
% - Includes a conditional to handle variables with uncertainties of 0. This is done by checking if an undefined value exists 
%   (resulting from a division by 0 in pder) and removing any variable associated with it. Variable are removed after the final value 
%   has already been calculated, so the answer is not affected. To learn more how this is done refer to the help function or Matlab docs.

try
    [value,fx] = pder(f,x,dx); % call pder, return final asnwer and partial derivatives
    nanArray = isnan(fx); 
    if ismember(true,nanArray) % if an undefined partial derivative exists
        indexOfNaN = find(nanArray); 
        fx(indexOfNaN) = []; 
        dx(indexOfNaN) = [];
        warning('''%s'' found to be errorless. Variable(s) ignored for uncertainty calculations\n',variableNames{indexOfNaN});
        variableNames(indexOfNaN) = []; 
        n = n - length(indexOfNaN); % decrease # of variables to represent only variables with uncertainties
    end
catch exception
    error('Possible reason(s):\n- function ''%s'' and variable(s) ''%s'' do not agree. \n- Problem with removing errorless variable.\n\nPlease try again.',functionString,variablesString);
end
totaldx = zeros(1,n); 
for i = 1:n
    totaldx(i) = fx(i) .* dx(i); % array of (partial derivatives) * (corresponding uncertainty)
end
df = sqrt(sum(totaldx.^2)); % total uncertainty
contributions = (totaldx./df).^2 * 100; % uncertainty contributions of each variable
variableLabels = categorical(variableNames,variableNames,'Ordinal',true); % x-axis labels
bar(variableLabels,contributions)
set(gca,'FontSize',16)
title('Sensitivity Analysis')
ylabel('Contribution to error [%]')
ylim([0 105])
for i = 1:n
    text(i,contributions(i),sprintf('%.2f %%',contributions(i)),'VerticalAlignment','bottom','HorizontalAlignment','center','Color','red','FontSize',14) % add value to each bar
end
fprintf('\n%s = %E +/- %E',desiredVariable,value,df) % print final answer


%% pder function
% - Calculates final value and partial derivatives

function [value,d] = pder(f,x,h)
    n = length(x); 
    switch n % switch case based on how many input variable
        case 1
            value = f(x(1));
            d(1) = (1/3) * ((4*(f(x(1)+h(1)) - f(x(1)-h(1)))./(2.*h(1))) - ((f(x(1)+(2*h(1))) - f(x(1)-(2*h(1))))./(4.*h(1))));
        case 2
            value = f(x(1),x(2));
            d(1) = (1/3) * ((4*(f(x(1)+h(1),x(2)) - f(x(1)-h(1),x(2)))./(2.*h(1))) - ((f(x(1)+(2*h(1)),x(2)) - f(x(1)-(2*h(1)),x(2)))./(4.*h(1))));
            d(2) = (1/3) * ((4*(f(x(1),x(2)+h(2)) - f(x(1),x(2)-h(2)))./(2.*h(2))) - ((f(x(1),x(2)+(2*h(2))) - f(x(1),x(2)-(2*h(2))))./(4.*h(2))));
        case 3
            value = f(x(1),x(2),x(3));
            d(1) = (1/3) * ((4*(f(x(1)+h(1),x(2),x(3)) - f(x(1)-h(1),x(2),x(3)))./(2.*h(1))) - ((f(x(1)+(2*h(1)),x(2),x(3)) - f(x(1)-(2*h(1)),x(2),x(3)))./(4.*h(1))));
            d(2) = (1/3) * ((4*(f(x(1),x(2)+h(2),x(3)) - f(x(1),x(2)-h(2),x(3)))./(2.*h(2))) - ((f(x(1),x(2)+(2*h(2)),x(3)) - f(x(1),x(2)-(2*h(2)),x(3)))./(4.*h(2))));
            d(3) = (1/3) * ((4*(f(x(1),x(2),x(3)+h(3)) - f(x(1),x(2),x(3)-h(3)))./(2.*h(3))) - ((f(x(1),x(2),x(3)+(2*h(3))) - f(x(1),x(2),x(3)-(2*h(3))))./(4.*h(3))));
        case 4
            value = f(x(1),x(2),x(3),x(4));
            d(1) = (1/3) * ((4*(f(x(1)+h(1),x(2),x(3),x(4)) - f(x(1)-h(1),x(2),x(3),x(4)))./(2.*h(1))) - ((f(x(1)+(2*h(1)),x(2),x(3),x(4)) - f(x(1)-(2*h(1)),x(2),x(3),x(4)))./(4.*h(1))));
            d(2) = (1/3) * ((4*(f(x(1),x(2)+h(2),x(3),x(4)) - f(x(1),x(2)-h(2),x(3),x(4)))./(2.*h(2))) - ((f(x(1),x(2)+(2*h(2)),x(3),x(4)) - f(x(1),x(2)-(2*h(2)),x(3),x(4)))./(4.*h(2))));
            d(3) = (1/3) * ((4*(f(x(1),x(2),x(3)+h(3),x(4)) - f(x(1),x(2),x(3)-h(3),x(4)))./(2.*h(3))) - ((f(x(1),x(2),x(3)+(2*h(3)),x(4)) - f(x(1),x(2),x(3)-(2*h(3)),x(4)))./(4.*h(3))));
            d(4) = (1/3) * ((4*(f(x(1),x(2),x(3),x(4)+h(4)) - f(x(1),x(2),x(3),x(4)-h(4)))./(2.*h(4))) - ((f(x(1),x(2),x(3),x(4)+(2*h(4))) - f(x(1),x(2),x(3),x(4)-(2*h(4))))./(4.*h(4))));
        case 5
            value = f(x(1),x(2),x(3),x(4),x(5));
            d(1) = (1/3) * ((4*(f(x(1)+h(1),x(2),x(3),x(4),x(5)) - f(x(1)-h(1),x(2),x(3),x(4),x(5)))./(2.*h(1))) - ((f(x(1)+(2*h(1)),x(2),x(3),x(4),x(5)) - f(x(1)-(2*h(1)),x(2),x(3),x(4),x(5)))./(4.*h(1))));
            d(2) = (1/3) * ((4*(f(x(1),x(2)+h(2),x(3),x(4),x(5)) - f(x(1),x(2)-h(2),x(3),x(4),x(5)))./(2.*h(2))) - ((f(x(1),x(2)+(2*h(2)),x(3),x(4),x(5)) - f(x(1),x(2)-(2*h(2)),x(3),x(4),x(5)))./(4.*h(2))));
            d(3) = (1/3) * ((4*(f(x(1),x(2),x(3)+h(3),x(4),x(5)) - f(x(1),x(2),x(3)-h(3),x(4),x(5)))./(2.*h(3))) - ((f(x(1),x(2),x(3)+(2*h(3)),x(4),x(5)) - f(x(1),x(2),x(3)-(2*h(3)),x(4),x(5)))./(4.*h(3))));
            d(4) = (1/3) * ((4*(f(x(1),x(2),x(3),x(4)+h(4),x(5)) - f(x(1),x(2),x(3),x(4)-h(4),x(5)))./(2.*h(4))) - ((f(x(1),x(2),x(3),x(4)+(2*h(4)),x(5)) - f(x(1),x(2),x(3),x(4)-(2*h(4)),x(5)))./(4.*h(4))));
            d(5) = (1/3) * ((4*(f(x(1),x(2),x(3),x(4),x(5)+h(5)) - f(x(1),x(2),x(3),x(4),x(5)-h(5)))./(2.*h(5))) - ((f(x(1),x(2),x(3),x(4),x(5)+(2*h(5))) - f(x(1),x(2),x(3),x(4),x(5)-(2*h(5))))./(4.*h(5))));
    end
end

%% getFunctionFromUser function
% - Prompts user to enter a function and extracts the required variables

function [f, functionString, desiredVariable, variableNames, n] = getFunctionFromUser()

    functionString = input('Enter your function in the format y=a*b\n', 's'); % prompt user to enter function
    equalSignPosition = strfind(functionString, '=');
    
    if length(equalSignPosition) ~= 1
        error('Possible reason(s):\n - Incorrect syntax used for your function ''%s'', there should be exactly one equal sign.', functionString)
    end
    
    formulaString = functionString(equalSignPosition+1:end); % extracts the formula from the user input (everything after the equal sign)
    
    % create array specifying which positions in user input contain letters, 
    % indicating a variable or a part of a variable (represented by a 1) 
    % and which positions do not, indicating a space or an operator (represented by a 0)
    variablePositions = isletter(functionString); 
    allVariables = strings(1, length(variablePositions)); % initialize empty array for variable names
    i = 1; % character counter
    j = 1; % variable counter
    
    % loop throuh variablePositions array and fill allVariables array
    while i <= length(variablePositions)
        if variablePositions(i) == 1 % if current indicie (i) represents a variable
            
            % loop through variablePositions array until a 0 appears (multi-character variable)
            while i <= length(variablePositions) && variablePositions(i) == 1
                allVariables(j) = strcat(allVariables(j), functionString(i)); % append subsequent characters to variable name
                i=i+1;
            end
        else % if current index (i) represents an operator, advance to next element in array
            i=i+1;
        end
        j=j+1;
    end

    allVariables = allVariables(~cellfun('isempty',allVariables)); % remove empty array elements
    
    desiredVariable = allVariables(1); % first element of allVariables array (variable being calculated)
    variableNames = allVariables(2:end); % all other elements of allVariables array (input variables)
    n = length(variableNames); % number of input variables
    
    variablesString = strjoin(cellstr(variableNames), ',');  % combine variables from array into one string
    
    try
        f = str2func([strcat('@(',variablesString,')') formulaString]); % convert variables and formula into function handle
    catch exception
        error('Possible reason(s):\n - Incorrect syntax used for your function ''%s''\n\nPlease try again.',functionString);
    end

end
