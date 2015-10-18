function newFun = expandFun(oldFun)

% get names of subfunctions and variables
funStr = func2str(oldFun);
funWorkspace = functions(oldFun);
funWorkspace = funWorkspace.workspace{1};
names = fieldnames(funWorkspace);

for i = 1:length(names)
    name = names{i};
    subWorkspace = funWorkspace.(name);
    lenName = length(name);
    lenFunStr = length(funStr);   
    c=1;
    
    % recursive loop for subfunctions
    if isa(subWorkspace,'function_handle')
        subFun = expandFun(subWorkspace);
        subStr = func2str(subFun);
        variablesLength = find(subStr==')');
        variablesLength = variablesLength(1) - 1;
        subStr = subStr(variablesLength+2:end);
        
        while c<lenFunStr-lenName+2
            
            % check name matches and is sorrounded by operators {+,-,*,/,^,(,)}
            if strcmp(name,funStr(c:c+lenName-1)) && opercmp(funStr,lenName,c)
                funStr = [funStr(1:c-1), '(', subStr, ')', funStr(c+lenName+variablesLength:end)];
                lenFunStr = length(funStr);
                c = c + length(subStr)+2;
            else
            c = c+1;
            end
        end
        
    % scalars and arrays    
    elseif isnumeric(subWorkspace)
        
        % scalars
        if length(subWorkspace)==1
            subStr = num2str(subWorkspace);
            
            while c<lenFunStr-lenName+2
                
                 % check name matches and is sorrounded by operators {+,-,*,/,^,(,)}
                if strcmp(name,funStr(c:c+lenName-1)) && opercmp(funStr,lenName,c)
                    funStr = [funStr(1:c-1), '(', subStr, ')', funStr(c+length(name):end)];
                    lenFunStr = length(funStr);
                    c = c + length(subStr)+2;
                else
                c = c+1;
                end
            end
            
        else
            
            while c<lenFunStr-lenName+2
                
                % check name matches and is sorrounded by operators {+,-,*,/,^,(,)}
                if strcmp(name,funStr(c:c+lenName-1)) && opercmp(funStr,lenName,c)
                    indLength = find(funStr(c+length(name):end)==')');
                    indLength = indLength(1)-2;
                    ind = str2num(funStr(c+lenName+1:c+lenName+indLength));
                    subStr = num2str(subWorkspace(ind));
                    funStr = [funStr(1:c-1), '(', subStr, ')', funStr(c+length(name)+indLength+2:end)];
                    lenFunStr = length(funStr);
                    c = c + length(subStr)+2;
                else
                c = c+1;
                end
            end
        end 
    end
end

newFun = str2func(funStr);

end

function bool = opercmp(funStr,lenName,c)

OperStr = {'+' ; '-' ; '*' ; '/' ; '^' ; '(' ; ')' ; '=' ; '<' ; '>'};

bool = 0;
if c <= 1
    suf = funStr(c+lenName);
    for s=1:7
        bool = bool + strcmp(OperStr{s},suf);
    end
    
elseif c+lenName > length(funStr)
    pre = funStr(c-1);
    for s=1:7
        bool = bool + strcmp(OperStr{s},pre);
    end
    
else
    pre = funStr(c-1);
    suf = funStr(c+lenName);
    for s=1:length(OperStr)
        bool = bool + strcmp(OperStr{s},pre) + strcmp(OperStr{s},suf);
    end
    bool = floor(bool/2);
end

end






