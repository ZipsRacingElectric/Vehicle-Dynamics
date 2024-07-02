% Save backup of model with <oldname><postfix_backup>
if (exist([Mdl_BaseName, postfix_backup, Mdl_FileExt]) ~= 4)
    [status,message,messageId] = copyfile (Mdl_FileName, [Mdl_BaseName, postfix_backup, Mdl_FileExt]);
    fprintf('Backuped: %s -> %s\n', Mdl_FileName, [Mdl_BaseName, postfix_backup, Mdl_FileExt]);
end