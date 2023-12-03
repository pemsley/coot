import { useState } from 'react'
// import './_App.css'
import './index.css';
import * as d3 from "d3";

function ToolButton(props) {
  return (
    <div className="button tool_button" onClick={props.onClick}>
      {props.caption}
      {/* {props.icon && <img src={props.icon} width="24px" /><br/>} */}
    </div>
  )
}

export function LhasaComponent() {
  const [lh, setLh] = useState(() => {
    // todo: figure out how to use this
    // return new Module.LhasaCanvas();
    return 0;
  });

  // function switch_tool(tool) {

  // };
  return (
    <>
      <div id="lhasa_editor">
        <div id="lhasa_hello_outer" className="horizontal_container">
          <img src="/icons/icons/hicolor_apps_scalable_coot-layla.svg" />
          <div id="lhasa_hello">
            <h3>Welcome to Lhasa!</h3>
            <p>
              Lhasa is a WebAssemby port of Layla - Coot's Ligand Editor.<br/>
              Lhasa is experimental software.
            </p>
            <p>
              This is a demo UI for development purposes.
            </p>
          </div>
        </div>
        <div id="molecule_tools_toolbar" className="horizontal_toolbar toolbar horizontal_container">
          {/* <div className="button tool_button" onclick="javascript:switch_tool(this,new Module.LhasaTransformTool(Module.LhasaTransformMode.Translation))" >Move</div>
          <div className="button tool_button" onclick="javascript:switch_tool(this,new Module.LhasaTransformTool(Module.LhasaTransformMode.Rotation))" >Rotate</div>
          <div className="button tool_button" onclick="javascript:switch_tool(this,new Module.LhasaFlipTool(Module.LhasaFlipMode.Horizontal))" >Flip around X</div>
          <div className="button tool_button" onclick="javascript:switch_tool(this,new Module.LhasaFlipTool(Module.LhasaFlipMode.Vertical))" >Flip around Y</div>
          <div className="button tool_button" onclick="javascript:switch_tool(this,new Module.LhasaRemoveHydrogensTool())" >Delete hydrogens</div>
          <div className="button tool_button" onclick="javascript:switch_tool(this,new Module.LhasaFormatTool())" >Format</div> */}
        </div>
        <div id="main_tools_toolbar" className="horizontal_toolbar toolbar horizontal_container">
          {/* <div className="button tool_button" onclick="javascript:switch_tool(this,new Module.LhasaBondModifier(Module.LhasaBondModifierMode.Single))" >Single Bond</div>
          <div className="button tool_button" onclick="javascript:switch_tool(this,new Module.LhasaBondModifier(Module.LhasaBondModifierMode.Double))" >Double Bond</div>
          <div className="button tool_button" onclick="javascript:switch_tool(this,new Module.LhasaBondModifier(Module.LhasaBondModifierMode.Triple))" >Triple Bond</div>
          <div className="button tool_button" onclick="javascript:switch_tool(this,new Module.LhasaGeometryModifier())" >Geometry</div>
          <div className="button tool_button" onclick="javascript:switch_tool(this,new Module.LhasaChargeModifier())" ><!--img src="icons/layla_charge_tool.svg" width="24px" / ><br/-->Charge</div>
          <div className="button tool_button" onclick="javascript:switch_tool(this,new Module.LhasaDeleteTool())" >Delete</div> */}
        </div>
        <div id="structure_toolbar" className="horizontal_toolbar toolbar horizontal_container">
          {/* <div className="button tool_button" onclick="javascript:switch_tool(this,new Module.LhasaStructureInsertion(Module.LhasaStructure.CycloPropaneRing))" >3-C</div>
          <div className="button tool_button" onclick="javascript:switch_tool(this,new Module.LhasaStructureInsertion(Module.LhasaStructure.CycloButaneRing))" >4-C</div>
          <div className="button tool_button" onclick="javascript:switch_tool(this,new Module.LhasaStructureInsertion(Module.LhasaStructure.CycloPentaneRing))" >5-C</div>
          <div className="button tool_button" onclick="javascript:switch_tool(this,new Module.LhasaStructureInsertion(Module.LhasaStructure.CycloHexaneRing))" >6-C</div>
          <div className="button tool_button" onclick="javascript:switch_tool(this,new Module.LhasaStructureInsertion(Module.LhasaStructure.BenzeneRing))" >6-Arom</div>
          <div className="button tool_button" onclick="javascript:switch_tool(this,new Module.LhasaStructureInsertion(Module.LhasaStructure.CycloHeptaneRing))" >7-C</div>
          <div className="button tool_button" onclick="javascript:switch_tool(this,new Module.LhasaStructureInsertion(Module.LhasaStructure.CycloOctaneRing))" >8-C</div> */}
        </div>
      </div>
    </>
  )
}

export function App() {
  // const [count, setCount] = useState(0)

  return (
    <>
     <LhasaComponent />
    </>
  )
}

// export default App
