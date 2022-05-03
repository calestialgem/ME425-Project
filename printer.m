classdef printer < handle
properties (Access = private)
    file
end
methods
function self = printer(file_name)
    [file, error_message] = fopen(file_name, 'w');
    if ~isempty(error_message)
        fprintf(2, "File %s Error Message: %s\n", file_name, error_message);
        return;
    end
    self.file = file;
end
function print(self, varargin)
    message = sprintf(varargin{:});
    fprintf(self.file, '%s\n', message);
    fprintf('%s\n', message);
end
function prmat(self, name, matrix, element)
    self.print("%s [%.0f, %.0f]: ", name, size(matrix, 1), size(matrix, 2));
    for k = 1 : size(matrix, 1)
        for j = 1 : size(matrix, 2)
            self.put(element, matrix(k, j));
        end
        self.put("\n");
    end
end
end
methods (Access = private)
function put(self, varargin)
    message = sprintf(varargin{:});
    fprintf(self.file, message);
    fprintf(message);
end
function delete(self)
    if fclose(self.file) ~= 0
        fprintf(2, "Error while closing the file\n");
        return;
    end
end
end
end
