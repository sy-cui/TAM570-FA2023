function savefig_pdf(filename);

    pos = get(gcf, 'Position');
    set(gcf,'PaperUnits', 'inches','PaperSize', [pos(3), pos(4)]);
    print(gcf, filename, '-dpdf');