// Minimal DOM stub sufficient to run app.js and observe what it does.
var CALLS = { vibrate: null, addAtoms: null, setStyle: null, animate: null };
function El(id, tag) {
  this.id=id; this.tagName=(tag||"DIV").toUpperCase(); this.children=[];
  this.classList={_s:{},add:function(c){this._s[c]=1},remove:function(c){delete this._s[c]},
                  toggle:function(c,v){ v?this._s[c]=1:delete this._s[c] },
                  contains:function(c){return !!this._s[c]}};
  this.dataset={}; this.style={}; this.checked=true; this.value=""; this.hidden=false;
  this.offsetTop=0; this.offsetHeight=18; this.scrollTop=0; this.clientHeight=190;
  var me=this;   // NOT `self`: dom.js defines a global `self`, which silently
                 // shadows any assumption that it means "this element"
  this.getBoundingClientRect=function(){
    if(me.id==="bondlist") return {top:400, bottom:400+me.clientHeight, height:me.clientHeight};
    var host=me._host;
    if(host) { var t = 400 - host.scrollTop + (me._rowTop||0);
               return {top:t, bottom:t+me.offsetHeight, height:me.offsetHeight}; }
    return {top:0, bottom:0, height:0};
  };
  this.textContent=""; this._listeners={};
  // innerHTML = "" must actually empty the node, as it does in a real DOM --
  // otherwise render()'s clear is a no-op and row counts accumulate.
  this._html="";
  Object.defineProperty(this,"innerHTML",{
    get:function(){return this._html;},
    set:function(v){ this._html=v; if(v==="") this.children.length=0; }
  });
  this.addEventListener=function(k,f){(this._listeners[k]=this._listeners[k]||[]).push(f);};
  this.fire=function(k,e){ (this._listeners[k]||[]).forEach(function(f){f(e||{});}); };
  this.appendChild=function(c){this.children.push(c); return c};
  this._qcache={};
  this.querySelector=function(sel){
    if(/pick/.test(sel)){ if(!this._pick){ this._pick=new El("input"); this._pick.className="pick";
      this._pick.dataset={i:String(this.dataset.index)}; } return this._pick; }
    // stable per selector: a fresh node each call means a listener attached
    // here is never the node a later lookup returns
    var m = /\.brow\[data-b="(\d+)"\]/.exec(sel);
    if(m && this.id==="bondlist"){
      // match on data-b, NOT on position: once the list is sorted the two
      // differ, and indexing by position would silently return another bond
      var rows=this.querySelectorAll(".brow"), want=m[1];
      for(var q=0;q<rows.length;q++) if(rows[q].dataset.b===want) return rows[q];
      return null;
    }
    if(!this._qcache[sel]) this._qcache[sel]=new El("q");
    return this._qcache[sel]; };
  this.querySelectorAll=function(sel){
    if(/brow/.test(sel) && (this.id==="bondlist" || this.id==="cmp-bonds")){
      // Cache on the rendered html: rebuilding the stubs on every call means a
      // listener attached here is never on the node a later lookup returns.
      if(this._rowsHtml !== this._html){
        var n=(this._html.match(/data-b="(\d+)"/g)||[]).map(function(m){
          return parseInt(m.replace(/\D/g,''),10); });
        // lay the rows out as a browser would, so scroll maths is exercised
        this._rows = n.map(function(i,k){ var e=new El("div"); e.className="brow";
          e.dataset={b:String(i)}; e.offsetHeight=18;
          // offsetTop as the DOM really reports it here: measured from the
          // positioned .viewer ancestor, not from the scrolling list
          e.offsetTop=900+k*18; e._rowTop=k*18; e._host=this; return e; }, this);
        this._rowsHtml = this._html;
      }
      return this._rows;
    }
    if(/th\.bsort/.test(sel) && this.id==="cmp-bonds"){
      if(this._thHtml !== this._html){
        var re=/<th class="bsort" data-bsort="([a-z]+)"(?: data-mode="(\d+)")?/g, mm, out=[];
        while((mm=re.exec(this._html))){
          var e=new El("th","TH"); e.className="bsort";
          e.dataset={bsort:mm[1]}; if(mm[2]!==undefined) e.dataset.mode=mm[2];
          out.push(e);
        }
        this._ths=out; this._thHtml=this._html;
      }
      return this._ths;
    }
    return []; };
  this.scrollIntoView=function(){};
  this.remove=function(){};
  this.click=function(){ this.fire('click',{target:{classList:{contains:function(){return false}}}}); };
}
var REG={};
function get(id){ if(!REG[id]) REG[id]=new El(id); return REG[id]; }
REG["payload"]=new El("payload"); REG["payload"].textContent=PAYLOAD_TXT;
REG["downloads"]=new El("downloads"); REG["downloads"].textContent=DL_TXT;
var TBODY=new El("tbody","TBODY");
var THEAD=["name","freq","Tx","Ty","Tz","Rx","Ry","Rz","V_S","label","irrep"].map(function(k){
  var th=new El("th","TH"); th.dataset={sort:k}; return th;
});
// The side panel's sort header, as result.html declares it.
var BHEAD=["pair","val"].map(function(k){
  var e=new El("bh","SPAN"); e.className="bsort"; e.dataset={bsort:k}; return e;
});
var document={
  getElementById:get,
  createElement:function(t){return new El("new",t)},
  querySelector:function(s){
    if (s.indexOf("tbody")>-1) return TBODY;
    var m=/\.bhead \.bsort\[data-bsort="([a-z]+)"\]/.exec(s);
    if (m) return BHEAD.filter(function(e){return e.dataset.bsort===m[1]})[0]||null;
    return new El("s"); },
  querySelectorAll:function(s){
    if (s.indexOf("tbody tr")>-1) return TBODY.children;
    if (s.indexOf("thead th")>-1) return THEAD;
    if (s===".bhead .bsort") return BHEAD;
    // markBondHeaders() sweeps every sortable header on the page at once --
    // the side panel's two and whatever the compare table currently shows
    if (s===".bsort") return BHEAD.concat(get("cmp-bonds").querySelectorAll("th.bsort"));
    // markBondRows() moves the compare table's selected-bond class without
    // re-rendering it, so it reaches the rows from the document down
    if (s==="#cmp-bonds tr.cbrow") return get("cmp-bonds").querySelectorAll("tr.cbrow");
    return [];
  },
  addEventListener:function(){},
  body:new El("body")
};
var window={matchMedia:function(){return{matches:false}}};
var URL={createObjectURL:function(){return "blob:x"},revokeObjectURL:function(){}};
function Blob(a,b){ this.parts=a; this.type=b.type; }
function getComputedStyle(){ return {getPropertyValue:function(){return "white"}}; }
// defaults the real inputs would carry
get("f-lo").value="0"; get("f-hi").value="99999";
get("amp").value="1"; get("frm").value="10";
get("f-ref").checked=true;
get("bonds").checked=true;
get("axes").checked=true;

["compare","cmp-n","cmp-scores","cmp-bonds","cmp-clear"].forEach(function(id){ get(id); });

var BONDROWS = true;
get("labels").checked=false;
