/* DOM stub rich enough for builder.js: real element tree, querySelectorAll
   with the class/:not selectors the builder actually uses. */
var window = this;
function El(tag) {
  this.tagName = (tag || "div").toUpperCase();
  this.children = []; this.className = ""; this._text = ""; this.value = "";
  this.type = ""; this.step = ""; this.placeholder = ""; this.title = "";
  this.hidden = false; this._l = {};
  var self = this;
  Object.defineProperty(this, "textContent", {
    get: function(){ return self._text; }, set: function(v){ self._text = v; }
  });
  Object.defineProperty(this, "innerHTML", {
    get: function(){ return ""; }, set: function(v){ if (v==="") self.children.length = 0; }
  });
  this.classList = {
    add: function(c){ if((" "+self.className+" ").indexOf(" "+c+" ")<0) self.className=(self.className+" "+c).trim(); },
    remove: function(c){ self.className=self.className.split(/\s+/).filter(function(x){return x!==c;}).join(" "); },
    contains: function(c){ return (" "+self.className+" ").indexOf(" "+c+" ")>=0; }
  };
  this.appendChild = function(c){ self.children.push(c); c.parent = self; return c; };
  this.insertBefore = function(c, ref){ var i=self.children.indexOf(ref); self.children.splice(i<0?self.children.length:i,0,c); return c; };
  this.remove = function(){ if(self.parent){ var i=self.parent.children.indexOf(self); if(i>=0) self.parent.children.splice(i,1); } };
  this.addEventListener = function(k,f){ self._l[k]=f; };
  this.click = function(){ if(self._l.click) self._l.click(); };
  this._all = function(out){ self.children.forEach(function(c){ out.push(c); c._all(out); }); return out; };
  this.querySelectorAll = function(sel){ return match(self._all([]), sel); };
  this.querySelector = function(sel){ return this.querySelectorAll(sel)[0] || null; };
}
function match(nodes, sel) {
  return nodes.filter(function (n) {
    return sel.split(",").some(function (part) {
      part = part.trim();
      var neg = null, m = part.match(/:not\(\.([\w-]+)\)/);
      if (m) { neg = m[1]; part = part.replace(m[0], ""); }
      var ok = part.split(/\s+/).filter(Boolean).every(function (tok) {
        if (tok.charAt(0) === ".") return tok.slice(1).split(".").every(function(c){ return n.classList.contains(c); });
        if (tok.indexOf(".") > 0) {
          var bits = tok.split(".");
          return n.tagName === bits[0].toUpperCase()
                 && bits.slice(1).every(function(c){ return n.classList.contains(c); });
        }
        return n.tagName === tok.toUpperCase();
      });
      return ok && (!neg || !n.classList.contains(neg));
    });
  });
}
var REG = {};
var document = {
  createElement: function(t){ return new El(t); },
  getElementById: function(id){ return REG[id] || null; }
};
["builder","b-atoms","b-modes","b-build","b-geom","b-bonds","b-modes-box","pasted","mode","form","b-view","b-vstat"]
  .forEach(function(id){ REG[id] = new El(id === "form" ? "form" : "div"); });
REG["b-atoms"].value = "3"; REG["b-modes"].value = "1"; REG["mode"].value = "paste";

document.body = new El("body");
function getComputedStyle(){ return { getPropertyValue: function(){ return "white"; } }; }
var CALLS = { models: [], style: null, zooms: 0 };
var $3Dmol = { createViewer: function(){ return {
  removeAllModels: function(){ CALLS.models = []; },
  addModel: function(data, fmt, opts){
    var n = parseInt(data.split("\n")[0], 10);
    var atoms = [];
    data.split("\n").slice(2).filter(function(l){return l.trim();}).forEach(function(l, i){
      var p = l.trim().split(/\s+/);
      atoms.push({ elem: p[0], x: +p[1], y: +p[2], z: +p[3],
                   bonds: [], bondOrder: [], index: null });
    });
    var m = { _a: atoms, selectedAtoms: function(){ return atoms; },
              setStyle: function(sel, st){ CALLS.style = st; },
              opts: opts };
    CALLS.models.push(m); return m;
  },
  zoomTo: function(){ CALLS.zooms++; }, render: function(){}
};}};
var setTimeout = function(f){ f(); return 0; };   // run the debounce immediately
var clearTimeout = function(){};
