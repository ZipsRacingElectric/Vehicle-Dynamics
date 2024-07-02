% Save html report which contains all changes
myformat = 'html'; % 'pdf','doc',...
if (publish_report==1)
    comparison = visdiff(strcat(trg_mdl, postfix_backup,  Mdl_FileExt),...
        strcat(trg_mdl, postfix_updated, Mdl_FileExt));
    publish(comparison,'Name',strcat(trg_mdl,'_changes'),'Format',myformat);
end