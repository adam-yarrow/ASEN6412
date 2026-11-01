function PrintCodeToPDF(targetDir)
    folderPath = fullfile(pwd,targetDir);
    mFiles = dir(fullfile(folderPath, '*.m'));
    
    outputDir = fullfile(folderPath,'codePdfs');
    
    for i = 1:length(mFiles)
        filePath = fullfile(folderPath, mFiles(i).name);
        options.format = 'pdf';
        options.outputDir = outputDir; 
        options.evalCode = false;
        publish(filePath, options);
        disp(['Published ', mFiles(i).name, ' to PDF.']);
    end
end