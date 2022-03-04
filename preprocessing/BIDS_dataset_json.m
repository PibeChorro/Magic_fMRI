%% Template Matlab script to create an BIDS compatible dataset_description.json file
% This example lists all required and optional fields.
% When adding additional metadata please use CamelCase
%
% Writing json files relies on the JSONio library
% https://github.com/gllmflndn/JSONio
% Make sure it is in the matab/octave path
%
% DHermes, 2017
% modified RG 201809
%%

function json_file = BIDS_dataset_json(rawdir)

dataset_description_json_name = fullfile(rawdir, 'dataset_description.json');

%% General fields, shared with MRI BIDS and MEG BIDS:

%%  Required fields:

dd_json.Name = 'Magic Project'; % name of the dataset

dd_json.BIDSVersion = '1.4.1'; % The version of the BIDS standard that was used

% The interpretation of the dataset. MUST be one of "raw" or "derivative". 
% For backwards compatibility, the default value is "raw".
dd_json.DatasetType = 'raw'; 

%%  Recommended fields:

dd_json.License = ''; % what license is this dataset distributed under? The
% use of license name abbreviations is suggested for specifying a license.
% A list of common licenses with suggested abbreviations can be found in appendix III.

dd_json.Authors = {'Vincent Plikat', 'Dr. Pablo Grassi', 'Prof. Dr. Andreas Bartels'}; % List of individuals who contributed to the
% creation/curation of the dataset

dd_json.Acknowledgements = "Dr. Michael Bannert and Jasper Bischoffberger"; % who should be acknowledge in helping to collect the data

dd_json.HowToAcknowledge = ''; % Instructions how researchers using this
% dataset should acknowledge the original authors. This field can also be used
% to define a publication that should be cited in publications that use the
% dataset.

dd_json.Funding = {'', '', ''}; % sources of funding (grant numbers)

% List of ethics committee approvals of the research protocols and/or protocol identifiers.
dd_json.EthicsApprovals = []; 

% a list of references to
% publication that contain information on the dataset, or links.
% dd_json.ReferencesAndLinks = {'', '', ''}; 

dd_json.DatasetDOI = ''; % the Document Object Identifier of the dataset
% (not the corresponding paper).

%% Write JSON

% this just makes the json file look prettier
% when opened in a text editor
json_options.indent = '    '; 

jsonSaveDir = fileparts(dataset_description_json_name);
if ~isdir(jsonSaveDir)
    fprintf('Warning: directory to save json file does not exist: %s \n', jsonSaveDir);
end

try
    jsonwrite(dataset_description_json_name, dd_json, json_options);
catch
    warning('%s\n%s\n%s\n%s', ...
            'Writing the JSON file seems to have failed.', ...
            'Make sure that the following library is in the matlab/octave path:', ...
            'https://github.com/gllmflndn/JSONio');
end
end
