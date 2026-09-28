const pptxgen = require('pptxgenjs');
const pres = new pptxgen();
pres.defineLayout({ name: 'PORTRAIT', width: 8.0, height: 11.0 });
pres.layout = 'PORTRAIT';
pres.title = 'Digital thread, detailed process';
const s = pres.addSlide();
s.background = { color: 'FFFFFF' };
const F = 'Times New Roman';
// Black and white (28 Sep 2026): every box white with a black outline, except the
// checks and outputs, which are light grey so the eye can find them.
const K = { in: ['FFFFFF', '000000'], geo: ['FFFFFF', '000000'], th: ['FFFFFF', '000000'],
            st: ['FFFFFF', '000000'], chk: ['E6E6E6', '000000'], out: ['E6E6E6', '000000'] };
// columns
const CX = 2.75, CW = 2.5, LX = 0.35, LW = 2.0, RX = 5.65, RW = 2.0;
const RH = 0.52, RG = 0.24, Y0 = 0.35;
const rowY = r => Y0 + r * (RH + RG);
function rect(x, y, w, h, txt, c, bold, fs) {
  s.addShape(pres.shapes.RECTANGLE, { x, y, w, h, fill: { color: c[0] }, line: { color: c[1], width: 1.25 } });
  s.addText(txt, { x, y, w, h, align: 'center', valign: 'middle', fontFace: F, fontSize: fs || 10.5,
    bold: !!bold, color: '000000', margin: 0.03, isTextBox: true });
}
function seg(x0, y0, x1, y1, arrow, color, dash) {
  s.addShape(pres.shapes.LINE, { x: Math.min(x0, x1), y: Math.min(y0, y1), w: Math.abs(x1 - x0), h: Math.abs(y1 - y0),
    flipH: x1 < x0, flipV: y1 < y0,
    line: { color: color || '000000', width: 1.25, dashType: dash || 'solid', endArrowType: arrow ? 'triangle' : 'none' } });
}
function path(pts, color, dash) {
  for (let i = 0; i < pts.length - 1; i++) seg(pts[i][0], pts[i][1], pts[i + 1][0], pts[i + 1][1], i === pts.length - 2, color, dash);
}
const C = (r, t, c, b) => rect(CX, rowY(r), CW, RH, t, c, b);
const L = (r, t, c) => rect(LX, rowY(r), LW, RH, t, c, false, 10);
const R = (r, t, c) => rect(RX, rowY(r), RW, RH, t, c, false, 10);
const mid = r => rowY(r) + RH / 2;
const inL = r => seg(LX + LW, mid(r), CX, mid(r), true);        // left input into centre
const inR = r => seg(RX, mid(r), CX + CW, mid(r), true);        // right input into centre
const outL = r => seg(CX, mid(r), LX + LW, mid(r), true);
const outR = r => seg(CX + CW, mid(r), RX, mid(r), true);

// row 0: three definitions
L(0, 'Aircraft missions\nJet-A fuel load', K.in);
C(0, 'Requirements\nCS-25 loads, hold time, pressures', K.in, true);
R(0, 'Material\nAl 2219-T87 allowables', K.in);
// row 1
L(1, 'Energy (LHV) equivalence\nfill 0.9', K.in);
C(1, 'Tank geometry and insulation\nR, L, heat budget Q, MLI', K.geo, true);
R(1, 'No-vent hold 14.7 h\nheat budget Q', K.in);
// row 2
L(2, 'Wall thickness t', K.geo);
C(2, 'Parametric shell model\nmesh seed ½√(Rt)', K.geo, true);
// row 3
L(3, 'Load cases\nBaseline, LC1, LC3, LC7', K.in);
C(3, 'Free surface\nplane normal to ĝ', K.geo, true);
R(3, 'Fill ratios\n0.10 to 0.90', K.in);
// row 4
L(4, 'Films h_wet, h_dry\nflux q = Q/A', K.th);
C(4, 'Thermal model\nwetted and dry wall', K.th, true);
R(4, 'Thermal properties\nconductivity κ', K.in);
// row 5
C(5, 'Steady-state heat transfer\n(Abaqus)', K.th, true);
// row 6
L(6, 'Check: ΔT = q / h_dry', K.chk);
C(6, 'Thermal results\ntemperature, heat flux', K.th, true);
R(6, 'Thermal design\nΔT, waterline flux', K.out);
// row 7
C(7, 'Temperature input\nconstant through thickness', K.th, true);
// row 8
L(8, 'Pressure + hydrostatic head\nshell inertia', K.st);
C(8, 'Structural model\ncold state stress-free', K.st, true);
R(8, 'Mechanical properties\nE, ν, α', K.in);
// row 9
L(9, 'Three rings\naft, mid anchor, fwd', K.st);
C(9, 'Linear static analysis\n(Abaqus)', K.st, true);
// row 10
L(10, 'Far-field stress\nreactions, displacement', K.out);
C(10, 'Structural results', K.st, true);
R(10, 'Yield at 1.5 bar\nbuckling at 0.187 bar', K.chk);
// row 11: decision
const dY = rowY(11) - 0.04, dH = RH + 0.18;
s.addShape(pres.shapes.DIAMOND, { x: CX + 0.25, y: dY, w: CW - 0.5, h: dH, fill: { color: 'FFFFFF' }, line: { color: '000000', width: 1.25 } });
s.addText('Wall just\nadequate?', { x: CX + 0.25, y: dY, w: CW - 0.5, h: dH, align: 'center', valign: 'middle',
  fontFace: F, fontSize: 10, bold: true, color: '000000', margin: 0, isTextBox: true });
// row 12
L(12, 'Mesh convergence\n3 levels', K.st);
C(12, 'Production and stability\nmatrices, 120 cases', K.st, true);
R(12, 'Verification\nhand checks, benchmark', K.chk);
// row 13
C(13, 'Comparison across\nthe three aircraft', K.out, true);

// spine
for (let r = 0; r < 13; r++) {
  if (r === 10) { seg(CX + CW / 2, rowY(10) + RH, CX + CW / 2, dY, true); continue; }
  if (r === 11) { seg(CX + CW / 2, dY + dH, CX + CW / 2, rowY(12), true); continue; }
  seg(CX + CW / 2, rowY(r) + RH, CX + CW / 2, rowY(r + 1), true);
}
// side inputs into the centre
inL(1); inR(1); inL(2); inL(3); inR(3); inL(4); inR(4); inL(8); inR(8); inL(9); inL(12);
// outputs and checks
outR(6); outL(6); outL(10); outR(10); outR(12);
// top row side boxes feed the row below
seg(LX + LW / 2, rowY(0) + RH, LX + LW / 2, rowY(1), true);
seg(RX + RW / 2, rowY(0) + RH, RX + RW / 2, rowY(1), true);
// decision: No -> back to wall thickness (left side loop); Yes -> down
const xl = LX - 0.18;
path([[CX + 0.25, dY + dH / 2], [xl, dY + dH / 2], [xl, mid(2)], [LX, mid(2)]], '000000');
s.addText('No: new t', { x: CX - 1.05, y: dY + dH / 2 - 0.26, w: 1.0, h: 0.22, fontFace: F, fontSize: 9.5, italic: true,
  color: '000000', margin: 0, align: 'right', isTextBox: true });
s.addText('Yes', { x: CX + CW / 2 + 0.12, y: dY + dH - 0.06, w: 0.5, h: 0.2, fontFace: F, fontSize: 9.5, italic: true,
  color: '000000', margin: 0, isTextBox: true });
pres.writeFile({ fileName: (process.argv[2] || 'DIAGRAMS_TO_REDRAW/08_digital_thread_detailed.pptx') }).then(() => console.log('ok'));
