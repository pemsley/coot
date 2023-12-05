import React from 'react'
import ReactDOM from 'react-dom/client'
import { App } from './Lhasa.tsx'
import './index.css'
import Module from '/lhasa.js?url'
import { MainModule } from './lhasa'

const Lhasa : MainModule = await Module();
export { Lhasa };

ReactDOM.createRoot(document.getElementById('root')!).render(
  <React.StrictMode>
    <App />
  </React.StrictMode>,
)
