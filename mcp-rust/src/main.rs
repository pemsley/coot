use serde_json::{json, Value};
use std::env;
use std::io::{Read, Write};
use std::net::TcpStream;
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::Duration;

const COOT_HOST: &str = "127.0.0.1"; // host stays constant; port is provided at runtime

// Atomic request id (mirrors itertools.count in the Python version)
static REQUEST_ID: AtomicU64 = AtomicU64::new(0);

fn send_coot_rpc(method: &str, params: Value, port: u16) -> Result<Value, String> {
    // Build JSON-RPC payload
    let id = REQUEST_ID.fetch_add(1, Ordering::SeqCst);
    let payload = json!({
        "jsonrpc": "2.0",
        "method": method,
        "params": params,
        "id": id,
    });

    let encoded = serde_json::to_vec(&payload).map_err(|e| format!("json serialize: {}", e))?;
    let len_prefix = (encoded.len() as u32).to_be_bytes();

    let addr = format!("{}:{}", COOT_HOST, port);
    let mut stream = TcpStream::connect(&addr).map_err(|e| format!("connect error: {}", e))?;

    // Optional: short timeouts to avoid blocking forever
    stream
        .set_read_timeout(Some(Duration::from_secs(10)))
        .ok();
    stream
        .set_write_timeout(Some(Duration::from_secs(10)))
        .ok();

    // Send length prefix + payload
    stream
        .write_all(&len_prefix)
        .map_err(|e| format!("write header error: {}", e))?;
    stream
        .write_all(&encoded)
        .map_err(|e| format!("write body error: {}", e))?;

    // Read 4-byte length prefix
    let mut header = [0u8; 4];
    stream
        .read_exact(&mut header)
        .map_err(|e| format!("read header error: {}", e))?;
    let body_len = u32::from_be_bytes(header) as usize;

    // Read exact body bytes
    let mut body_buf = vec![0u8; body_len];
    stream
        .read_exact(&mut body_buf)
        .map_err(|e| format!("read body error: {}", e))?;

    let response: Value = serde_json::from_slice(&body_buf).map_err(|e| format!("json parse: {}", e))?;
    Ok(response)
}

/// Send Python code to Coot using the "python.exec" method.
///
/// The `code` string can be multi-line Python source; it will be placed
/// in the JSON params as: { "code": "<your python source>" }.
fn run_python(code: &str, port: u16) -> Result<Value, String> {
    let params = json!({ "code": code });
    send_coot_rpc("python.exec", params, port)
}

/// Convenience: list available tools
fn list_available_tools(port: u16) -> Result<Value, String> {
    let params = json!({});
    send_coot_rpc("mcp.list_tools", params, port)
}

fn print_usage_and_exit(program: &str) -> ! {
    eprintln!("Usage: {} [PORT]", program);
    eprintln!("If PORT is omitted, defaults to 9090.");
    std::process::exit(1);
}

fn main() {
    // Parse optional positional port argument
    let mut args = env::args();
    let program = args.next().unwrap_or_else(|| "coot_socket_bridge".to_string());
    let port: u16 = match args.next() {
        Some(s) => match s.parse::<u16>() {
            Ok(p) => p,
            Err(_) => {
                eprintln!("Invalid port number: {}", s);
                print_usage_and_exit(&program);
            }
        },
        None => 9090,
    };

    println!("Using Coot at {}:{}", COOT_HOST, port);

    // Example: send a short python expression
    let python_code = r#"
# simple example
result = 1 + 2
result
"#;

    match run_python(python_code, port) {
        Ok(resp) => println!("run_python response: {}", resp),
        Err(e) => eprintln!("run_python error: {}", e),
    }

    // Example: list tools
    match list_available_tools(port) {
        Ok(resp) => println!("list_available_tools response: {}", resp),
        Err(e) => eprintln!("list_available_tools error: {}", e),
    }
}