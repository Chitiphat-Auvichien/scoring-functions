/* Shadow the $3Dmol global with a proxy whose createViewer returns a
 * headless viewer. Model creation still goes through the REAL GLModel, so
 * parsing, bonding and vibrate() are genuinely exercised -- only WebGL
 * rendering is stubbed. */
var CALLS = { vibrate:null, animate:null, style:null, viewer:null,
              animateCount:0, stopCount:0, live:0, zooms:0, frame:null,
              arrows:[], labels:[], cylinders:[] };
var REAL = $3Dmol;
$3Dmol = {
  GLModel: REAL.GLModel,
  Parsers: REAL.Parsers,
  createViewer: function () {
    var v = {
      models: [],
      getModelOpt: function (o) { return o || {}; },
      removeAllModels: function () { this.models = []; },
      removeAllShapes: function () { CALLS.cylinders = []; CALLS.arrows = []; },
      addModel: function (data, fmt, opts) {
        var m = new REAL.GLModel(this.models.length, opts || {}, null);
        m.addMolData(data, fmt, opts || {});
        var _v = m.vibrate.bind(m), _s = m.setStyle.bind(m);
        m.vibrate = function (n, a, b, ar) {
          CALLS.vibrate = { frames:n, amp:a, both:b, arrow:ar }; return _v(n,a,b,ar);
        };
        m.setStyle = function (sel, st) { CALLS.style = st; return _s(sel, st); };
        this.models.push(m); return m;
      },
      addArrow: function (spec) { CALLS.arrows.push(spec); },
      addCylinder: function (spec) { CALLS.cylinders.push(spec); },
      addLabel: function (t, spec) { CALLS.labels.push({text: t, spec: spec}); },
      removeAllLabels: function () { CALLS.labels = []; CALLS.arrows = []; },
      resize: function () { CALLS.resizes = (CALLS.resizes||0) + 1; },
      linkedViewers: [],
      linkViewer: function (o) { this.linkedViewers.push(o); return this; },
      getView: function () { return this._view || [0,0,0,1,0,0,0,1]; },
      setView: function (v) { CALLS.setViews = (CALLS.setViews||0)+1; this._view = v; },
      zoomTo: function () { CALLS.zooms++; }, render: function () {},
      // count loops started vs stopped, and model frame index
      animate: function (o) { CALLS.animate = o; CALLS.animateCount++; CALLS.live++; },
      stopAnimate: function () { CALLS.stopCount++; CALLS.live = 0; },
      setFrame: function (n) { CALLS.frame = n; return { then: function (f) { f(); } }; }
    };
    CALLS.viewer = v; (CALLS.viewers = CALLS.viewers || []).push(v); return v;
  }
};
