% Copyright (c) 2017-2018 Paul Irofti <paul@irofti.net>
% 
% Permission to use, copy, modify, and/or distribute this software for any
% purpose with or without fee is hereby granted, provided that the above
% copyright notice and this permission notice appear in all copies.
% 
% THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
% WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
% MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
% ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
% WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
% ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
% OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.

function mkpdf(figdir, figname)
    cfig = sprintf('%s%s', figdir, figname);
    set(gca(), 'LooseInset', get(gca(), 'TightInset'));
    f = gcf;
    set(f, 'PaperUnits','centimeters');
    set(f, 'Units','centimeters');
    pos=get(f,'Position');
    set(f, 'PaperSize', [pos(3) pos(4)]);
    set(f, 'PaperPositionMode', 'manual');
    set(f, 'PaperPosition',[0 0 pos(3) pos(4)]);
    hgexport(f, [cfig '.pdf'], ...
        hgexport('factorystyle'), 'Format', 'pdf');
end